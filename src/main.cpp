#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <sstream>
#include <thread>
#include <cassert>

#include <boost/any.hpp>
#include <boost/filesystem.hpp>

#include <pbbam/BamReader.h>
#include <pbbam/BamRecord.h>
#include <pbbam/BamWriter.h>

#include <Complete-Striped-Smith-Waterman-Library/src/ssw_cpp.h>

#include "common.hpp"
#include "version.inc"
#include "threads.hpp"

using namespace std;
using namespace PacBio::BAM;

enum Arguments {
    PRIMER
    , OUTPUT
    , THREADS
    , BULKSIZE
    , INPUT
    , MIN_SW_SCORE
    , MIN_SW_SCORE_DIFF
    , SW_MATCH_SCORE
    , SW_MISMATCH_PENALTY
    , SW_GAP_OPEN_PENALTY
    , SW_GAP_EXT_PENALTY
    , SIZE
};

using argument_type = array<string, Arguments::SIZE>;

template <bool is_left>
void SplitBam(const BamRecord& record
              , BamRecord& outbam
              , const StringView& run_name
              , const StringView& zmw
              , int left_start
              , int right_end
              , const StripedSmithWaterman::Alignment& alignment
             ) {
    // name
    stringstream name;
    name << run_name.to_string() << '/' << zmw << '/';
    if (is_left) {
        name << left_start << '_' << (left_start + alignment.ref_begin);
    } else {
        name << (left_start + alignment.ref_end + 1) << '_' << right_end;
    }
    outbam.Impl().Name(name.str());
    // sequence
    outbam.Impl().SetSequenceAndQualities(
        is_left ?
        record.Sequence().substr(0, alignment.ref_begin) :
        record.Sequence().substr(alignment.ref_end + 1)
        , ""
    );
    // tags
    const auto& original_tags = record.Impl().Tags();
    TagCollection newtags;
    // direct copy tags
    for (auto t  : {"RG", "np", "rq", "sn", "zm"}) {
        newtags[t] = original_tags.at(t);
    }
    // tags need to substr
    for (auto t : {"dq", "dt", "iq", "mq", "sq"}) {
        if (is_left)
            newtags[t] = original_tags.at(t).ToString().substr(0, alignment.ref_begin);
        else
            newtags[t] = original_tags.at(t).ToString().substr(alignment.ref_end + 1);
    }
    // special tag
    newtags["qs"] = (is_left ? left_start : left_start + alignment.ref_end + 1);
    newtags["qe"] = (is_left ? left_start + alignment.ref_begin : right_end);
    newtags["cx"] = original_tags.at("cx").ToUInt8()
        | (is_left ? PacBio::BAM::LocalContextFlags::ADAPTER_BEFORE : PacBio::BAM::LocalContextFlags::ADAPTER_AFTER);
    // TODO: this can be more efficient
    auto vec = record.IPD().Encode();
    if (is_left) {
        newtags["ip"] = decltype(vec){vec.cbegin(), vec.cbegin() + alignment.ref_begin};
    } else {
        newtags["ip"] = decltype(vec){vec.cbegin() + alignment.ref_end + 1, vec.cend()};
    }
    outbam.Impl().Tags(newtags);
}


std::mutex k_io_mx;

class BamSplitter {
private:
    using queue_type = MultiThreadSafeQueue<vector, BamRecord>;

    uint16_t min_sw_score_;
    uint16_t min_sw_diff_;
    uint8_t match_score_;
    uint8_t mismatch_penalty_;
    uint8_t gap_open_penalty_;
    uint8_t gap_ext_penalty_;
    queue_type& queue_;
    BamWriter& writer_;
    const BamHeader& header_;
    const string& primer_seq_;

public:
    BamSplitter(queue_type& q
                , BamWriter& w
                , const string& p
                , const BamHeader& h
                , uint16_t min_sw_score
                , uint16_t min_sw_diff
                , uint8_t match_score
                , uint8_t mismatch_penalty
                , uint8_t gap_open_penalty
                , uint8_t gap_ext_penalty
               )
        : queue_(q)
          , writer_(w)
          , primer_seq_(p)
          , header_(h)
          , min_sw_score_(min_sw_score)
          , min_sw_diff_(min_sw_diff)
          , match_score_(match_score)
          , mismatch_penalty_(mismatch_penalty)
          , gap_open_penalty_(gap_open_penalty)
          , gap_ext_penalty_(gap_ext_penalty) {}

    BamSplitter(const BamSplitter&) = delete;

    BamSplitter(BamSplitter&& other)
        :
        queue_(other.queue_)
        , writer_(other.writer_)
        , header_(other.header_)
        , primer_seq_(other.primer_seq_)
        , min_sw_score_(other.min_sw_score_)
        , min_sw_diff_(other.min_sw_diff_)
        , match_score_(other.match_score_)
        , mismatch_penalty_(other.mismatch_penalty_)
        , gap_open_penalty_(other.gap_open_penalty_)
        , gap_ext_penalty_(other.gap_ext_penalty_) {}

    BamSplitter& operator=(const BamSplitter&) = delete;

    void operator()() {
        vector<BamRecord> outputs;
        outputs.reserve(queue_.Capacity());
        string fullname;
        int left_start, right_end;
        StripedSmithWaterman::Aligner aligner{match_score_, mismatch_penalty_, gap_open_penalty_, gap_ext_penalty_};
        StripedSmithWaterman::Filter filter;
        StripedSmithWaterman::Alignment alignment;
        // begin process data
        auto data = queue_.FillAndPop();
        while (!data.empty()) {
            for (auto& record : data) {
                alignment.Clear();
                aligner.Align(primer_seq_.c_str()
                              , record.Sequence().c_str()
                              , record.Sequence().size()
                              , filter
                              , &alignment
                );
                // filter
                if (alignment.sw_score < min_sw_score_
                    || alignment.sw_score - alignment.sw_score_next_best < min_sw_diff_) {
                    #ifndef NDEBUG
                    if (alignment.sw_score < min_sw_score_) {
                        fprintf(stderr
                                , "[1]\t%d\t%d\t%s\n"
                                , alignment.sw_score
                                , alignment.sw_score_next_best
                                , record.Sequence().c_str());
                    } else {
                        fprintf(stderr
                                , "[2]\t%d\t%d\t%s\n"
                                , alignment.sw_score
                                , alignment.sw_score_next_best
                                , record.Sequence().c_str());
                    }
                    #endif
                    continue;
                }
                // fix name
                fullname = record.FullName();
                auto tokens = Utils::Tokenize(fullname, '/');
                auto tokens2 = Utils::Tokenize(tokens[2], '_');
                if (!Utils::StringViewTo(tokens2[0], left_start) || !Utils::StringViewTo(tokens2[1], right_end)) {
                    Utils::Error("failed to convert start or end");
                }
                // fix sequence
                if (alignment.ref_begin > 0) {
                    BamRecord l(header_);
                    SplitBam<true>(record, l, tokens[0], tokens[1], left_start, right_end, alignment);
                    outputs.push_back(std::move(l));
                }
                if (left_start + alignment.ref_end + 1 < right_end) {
                    BamRecord r(header_);
                    SplitBam<false>(record, r, tokens[0], tokens[1], left_start, right_end, alignment);
                    outputs.push_back(std::move(r));
                }
            } // end of processing each BamRecord from queue
            {
                lock_guard<mutex> lock(k_io_mx);
                for (const auto& o : outputs) {
                    writer_.Write(o);
                }
            }
            outputs.clear();
            data = queue_.FillAndPop();
        }
    }
};


int SplitterMT(const argument_type& args) {
    auto out_file_name = args[Arguments::OUTPUT];
    auto primer_seq = args[Arguments::PRIMER];
    auto min_sw_score = static_cast<uint16_t>(stoi(args[Arguments::MIN_SW_SCORE]));
    auto max_sw_diff = static_cast<uint16_t>(stoi(args[Arguments::MIN_SW_SCORE_DIFF]));
    auto match_score = static_cast<uint8_t>(stoi(args[Arguments::SW_MATCH_SCORE]));
    auto mismatch_penalty = static_cast<uint8_t>(stoi(args[Arguments::SW_MISMATCH_PENALTY]));
    auto gap_open_penalty = static_cast<uint8_t>(stoi(args[Arguments::SW_GAP_OPEN_PENALTY]));
    auto gap_ext_penalty = static_cast<uint8_t>(stoi(args[Arguments::SW_GAP_EXT_PENALTY]));
    auto subread_bam_file = args[Arguments::INPUT];
    BamReader subread_bam_fh(subread_bam_file);
    MultiThreadSafeQueue<vector, BamRecord> queue(subread_bam_fh, stoul(args[Arguments::BULKSIZE]));
    auto header = subread_bam_fh.Header().DeepCopy();
    BamWriter out_fh(out_file_name
                     , header
                     , BamWriter::CompressionLevel::CompressionLevel_4
                     , 1
                     , BamWriter::BinCalculation_OFF
    );

    int numThreads = stoi(args[Arguments::THREADS]);
    vector<thread> threads;
    for (int i = 0; i < numThreads; ++i) {
        threads.emplace_back(BamSplitter{queue, out_fh, primer_seq, header, min_sw_score, max_sw_diff, match_score
                                         , mismatch_penalty, gap_open_penalty, gap_ext_penalty});
    }
    for (auto& t : threads) {
        if (t.joinable()) {
            t.join();
        }
    }
    return EXIT_SUCCESS;
}

argument_type ArgumentParse(int argc, char** argv) {
    string usage =
        KERNAL_GREEN
            "This program split pacBio bam files based on a custom smrtbell sequence.\n"
            "version: v"
            PROGRAM_VERSION
            "\nusage: [options] input.bam\n\n"
            "options:\n"
            KERNAL_RED "\n[required]\n"
            "\tinput.bam\n"
            KERNAL_CYAN
            "\n[optional]\n"
            "\t-o      output bam filename, if not provided, the prefix of input bam + refarm.bam will be used\n"
            "\t-p      primer sequence, default: " DEFAULT_PRIMER_SEQ "\n"
            "\t-t      number of threads to use, default: " DEFAULT_NUM_THREADS "\n"
            KERNAL_YELLOW
            "\n[advanced]\n"
            "\t-b      bulk of records sent to each thread every time, default: " DEFAULT_BULK_SIZE "\n"
            "\t-m      minimal Smith-Waterman score between read and adaptor, default: " DEFAULT_MIN_SW_SCORE "\n"
            "\t-f      minimal Smith-Waterman score allowed between best and second-best alignments, default: " DEFAULT_MIN_SW_DIFF "\n"
            "\t-M      Score for a match, default: " DEFAULT_SW_MATCH_SCORE "\n"
            "\t-S      Penalty for a mismatch, default: " DEFAULT_SW_MISMATCH_PENALTY "\n"
            "\t-O      Penalty for a gap opening, default: " DEFAULT_SW_GAP_OPEN_PENALTY "\n"
            "\t-E      Penalty for a gap extension, default: " DEFAULT_SW_GAP_EXT_PENALTY "\n"
            KERNAL_RESET;

    argument_type arguments;
    int c;
    while ((c = getopt(argc, argv, "p:o:t:b:f:m:M:S:O:E:h")) != -1) {
        switch (c) {
            case 'p':arguments[Arguments::PRIMER] = optarg;
                break;
            case 'o':arguments[Arguments::OUTPUT] = optarg;
                break;
            case 't':arguments[Arguments::THREADS] = optarg;
                break;
            case 'b':arguments[Arguments::BULKSIZE] = optarg;
                break;
            case 'f':arguments[Arguments::MIN_SW_SCORE_DIFF] = optarg;
                break;
            case 'm':arguments[Arguments::MIN_SW_SCORE] = optarg;
                break;
            case 'M':arguments[Arguments::SW_MATCH_SCORE] = optarg;
                break;
            case 'S':arguments[Arguments::SW_MISMATCH_PENALTY] = optarg;
                break;
            case 'O':arguments[Arguments::SW_GAP_OPEN_PENALTY] = optarg;
                break;
            case 'E':arguments[Arguments::SW_GAP_EXT_PENALTY] = optarg;
                break;
            case 'h':
            default:cerr << usage;
                exit(EXIT_FAILURE);
        }
    }
    if (optind == argc) {
        Utils::Error("Please provide input bam files");
    }
    arguments[Arguments::INPUT] = argv[optind];
    if (access(arguments[Arguments::INPUT].c_str(), F_OK) == -1) {
        Utils::Error("Input file " + arguments[Arguments::INPUT] + " does not exist");
    }
    if (arguments[Arguments::OUTPUT].empty()) {
        string prefix = boost::filesystem::basename(arguments[Arguments::INPUT]);
        arguments[Arguments::OUTPUT] = prefix.substr(0, prefix.rfind(".bam")) + ".refarm.bam";
        Utils::Warning(
            "The user has not provided a output file prefix using -o option, will use the prefix of input bam"
        );
    }
    if (arguments[Arguments::PRIMER].empty()) { arguments[Arguments::PRIMER] = DEFAULT_PRIMER_SEQ; }
    if (arguments[Arguments::THREADS].empty()) { arguments[Arguments::THREADS] = DEFAULT_NUM_THREADS; }
    if (arguments[Arguments::BULKSIZE].empty()) { arguments[Arguments::BULKSIZE] = DEFAULT_BULK_SIZE; }
    if (arguments[Arguments::MIN_SW_SCORE].empty()) { arguments[Arguments::MIN_SW_SCORE] = DEFAULT_MIN_SW_SCORE; }
    if (arguments[Arguments::MIN_SW_SCORE_DIFF].empty()) {
        arguments[Arguments::MIN_SW_SCORE_DIFF] = DEFAULT_MIN_SW_DIFF;
    }
    if (arguments[Arguments::SW_MATCH_SCORE].empty()) { arguments[Arguments::SW_MATCH_SCORE] = DEFAULT_SW_MATCH_SCORE; }
    if (arguments[Arguments::SW_MISMATCH_PENALTY].empty()) {
        arguments[Arguments::SW_MISMATCH_PENALTY] = DEFAULT_SW_MISMATCH_PENALTY;
    }
    if (arguments[Arguments::SW_GAP_OPEN_PENALTY].empty()) {
        arguments[Arguments::SW_GAP_OPEN_PENALTY] = DEFAULT_SW_GAP_OPEN_PENALTY;
    }
    if (arguments[Arguments::SW_GAP_EXT_PENALTY].empty()) {
        arguments[Arguments::SW_GAP_EXT_PENALTY] = DEFAULT_SW_GAP_EXT_PENALTY;
    }
    return arguments;
}

int main(int argc, char** argv) {
    auto args = ArgumentParse(argc, argv);
    return SplitterMT(args);
}
