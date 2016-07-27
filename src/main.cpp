#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <sstream>
#include <thread>
#include <cstdlib>
#include <cassert>
#include <boost/any.hpp>
#include <boost/filesystem.hpp>

#include <pbbam/BamFile.h>
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
    , SIZE
};
using argument_type = array<string, Arguments::SIZE>;

#ifndef DEFAULT_PRIMER_SEQ
#define DEFAULT_PRIMER_SEQ "ATCTCTCTCAATTTTTTTTTTTTTTTTTTTTTTTAAGAGAGAGAT"
#endif

#ifndef DEFAULT_NO_THREADS
#define DEFAULT_NO_THREADS "8"
#endif

#ifndef DEFAULT_BULK_SIZE
#define DEFAULT_BULK_SIZE "100"
#endif


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
    }
    else {
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

    queue_type& queue_;
    BamWriter& writer_;
    const BamHeader& header_;
    const string& primer_seq_;
public:
    BamSplitter(queue_type& q, BamWriter& w, const string& p, const BamHeader& h)
      : queue_(q)
        , writer_(w)
        , primer_seq_(p)
        , header_(h) {}

    BamSplitter(const BamSplitter&) = delete;

    BamSplitter(BamSplitter&& other)
      :
      queue_(other.queue_)
      , writer_(other.writer_)
      , header_(other.header_)
      , primer_seq_(other.primer_seq_) {

    }

    BamSplitter& operator=(const BamSplitter&) = delete;

    void operator()() {
        vector<BamRecord> outputs;
        outputs.reserve(queue_.Capacity());
        string fullname;
        int left_start, right_end;
        StripedSmithWaterman::Aligner aligner;
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
    string out_file_name = args[Arguments::OUTPUT];
    string primer_seq = args[Arguments::PRIMER];
    auto subread_bam_file = args[Arguments::INPUT];
    BamReader subread_bam_fh(subread_bam_file);
    MultiThreadSafeQueue<vector, BamRecord> queue(subread_bam_fh, stoi(args[Arguments::BULKSIZE]));
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
        threads.emplace_back(BamSplitter{queue, out_fh, primer_seq, header});
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

        KERNAL_CYAN
        "\n[optional]\n"
        "\t-o      output bam filename, if not provided, the prefix of input bam + refarm.bam will be used\n"
        "\t-p      primer sequence, default: " DEFAULT_PRIMER_SEQ "\n"
        "\t-t      number of threads to use, default: " DEFAULT_NO_THREADS "\n"

        KERNAL_YELLOW
        "\n[advanced]\n"
        "\t-b      bulk of records sent to each thread every time, default: " DEFAULT_BULK_SIZE "\n"
        KERNAL_RESET;

    argument_type arguments;
    int c;
    while ((c = getopt(argc, argv, "p:o:t:b:h")) != -1) {
        switch (c) {
            case 'p':
                arguments[Arguments::PRIMER] = optarg;
                break;
            case 'o':
                arguments[Arguments::OUTPUT] = optarg;
                break;
            case 't':
                arguments[Arguments::THREADS] = optarg;
                break;
            case 'b':
                arguments[Arguments::BULKSIZE] = optarg;
                break;
            case 'h':
            default:
                cerr << usage;
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

    if (arguments[Arguments::PRIMER].empty()) {
        arguments[Arguments::PRIMER] = DEFAULT_PRIMER_SEQ;
    }
    if (arguments[Arguments::THREADS].empty()) {
        arguments[Arguments::THREADS] = DEFAULT_NO_THREADS;
    }
    if (arguments[Arguments::BULKSIZE].empty()) {
        arguments[Arguments::BULKSIZE] = DEFAULT_BULK_SIZE;
    }
    if (arguments[Arguments::OUTPUT].empty()) {
        string prefix = boost::filesystem::basename(arguments[Arguments::INPUT]);
        arguments[Arguments::OUTPUT] = prefix.substr(0, prefix.rfind(".bam")) + ".refarm.bam";
        Utils::Warning(
          "The user has not provided a output file prefix using -o option, will use the prefix of input bam"
                      );
    }
    return arguments;
}

int main(int argc, char** argv) {
    auto args = ArgumentParse(argc, argv);
    return SplitterMT(args);
}
