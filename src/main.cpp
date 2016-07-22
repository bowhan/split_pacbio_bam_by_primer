#include <iostream>
#include <string>
#include <sstream>
#include <cassert>
#include <pbcopper/cli/CLI.h>
#include <pbbam/BamFile.h>
#include <pbbam/BamReader.h>
#include <pbbam/BamRecord.h>
#include <pbbam/BamWriter.h>
#include <Complete-Striped-Smith-Waterman-Library/src/ssw_cpp.h>

#include "common.hpp"
#include "config.h"

using namespace std;

namespace {
using namespace PacBio::CLI;
using namespace PacBio::BAM;

static Interface CreateCLI() {
    Interface i{"remove_adaptor_from_subreads_bam", "a tool to remove custom smrtbell adaptor from subreads bam"
                , PROGRAM_VERSION};
    i.AddHelpOption().AddVersionOption().AddPositionalArgument({"input", "Input subread.bam file"});
    i.AddOptions({
                   {"output", {"o", "output"}, "ouput file prefix", Option::StringType("")}
                   , {"primer", {"p", "primers"}, "primer sequence in fasta format", Option::StringType("")}
                   , {"threads", {"t", "threads"}, "Threads used", Option::IntType(8)}
                 }
                );
    return i;
}

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
    // flags
//    outbam.Impl().Flag(record.Impl().Flag());
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
    newtags["ip"] = record.IPD().Encode();

    outbam.Impl().Tags(newtags);
}

/// \brief
/// \param args
/// \return
int Runner(const PacBio::CLI::Results& args) {
    string outputPrefix = args["output"];
    string primer_seq = args["primer"];

    if (outputPrefix.empty()) {
        cerr << "Error: please provide the output prefix" << endl;
        return EXIT_FAILURE;
    }
    auto subread_bam_files = args.PositionalArguments();

    int zmw, left_start, right_end;
    // align
    // Declares a default Aligner
    StripedSmithWaterman::Aligner aligner;
    // Declares a default filter
    StripedSmithWaterman::Filter filter;
    // Declares an alignment that stores the result
    StripedSmithWaterman::Alignment alignment;

    for (const auto& subread_bam_file : subread_bam_files) {
        BamReader subread_bam_fh(subread_bam_file);
        string out_file_name = subread_bam_file.substr(0, subread_bam_file.rfind(".bam")) + ".refarm.bam";
        auto header = subread_bam_fh.Header().DeepCopy();
        BamWriter out_fh(out_file_name
                         , header
                         , BamWriter::CompressionLevel::CompressionLevel_4
                         , args["threads"]
                         , BamWriter::BinCalculation_OFF
                        );
        BamRecord record;
        string fullname; // for BAM record
        while (subread_bam_fh.GetNext(record)) {
            // out
            BamRecord l(header), r(header);
            // Aligns the query to the ref
            alignment.Clear();
            aligner.Align(primer_seq.c_str(), record.Sequence().c_str(), record.Sequence().size(), filter, &alignment);
            // fix name
            fullname = record.FullName();
            auto tokens = Utils::Tokenize(fullname, '/');
            auto tokens2 = Utils::Tokenize(tokens[2], '_');
            if (!Utils::StringViewTo(tokens2[0], left_start) || !Utils::StringViewTo(tokens2[1], right_end)) {
                std::cerr << "failed to convert start or end\n";
                return 1;
            }
            // fix sequence
            SplitBam<true>(record, l, tokens[0], tokens[1], left_start, right_end, alignment);
            SplitBam<false>(record, r, tokens[0], tokens[1], left_start, right_end, alignment);
            // write
            out_fh.Write(l);
            out_fh.Write(r);
        }
    }
}

} // end of namespace

int main(int argc, char** argv) {
    return PacBio::CLI::Run(argc, argv, CreateCLI(), &Runner);
}