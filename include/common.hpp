#ifndef SPLIT_POLYT_FROM_PBBAM_COMMON_HPP
#define SPLIT_POLYT_FROM_PBBAM_COMMON_HPP

#include <vector>
#include <sstream>
#include <boost/utility/string_ref.hpp>
#include "kernel_color.hpp"

#ifndef DEFAULT_PRIMER_SEQ
#define DEFAULT_PRIMER_SEQ "ATCTCTCTCAATTTTTTTTTTTTTTTTTTTTTTTAAGAGAGAGAT"
#endif

#ifndef DEFAULT_NUM_THREADS
#define DEFAULT_NUM_THREADS "4"
#endif

#ifndef DEFAULT_BULK_SIZE
#define DEFAULT_BULK_SIZE "500"
#endif

#ifndef DEFAULT_MIN_LEN_REPORT
#define DEFAULT_MIN_LEN_REPORT "100"
#endif

#ifndef DEFAULT_MIN_SW_SCORE
#define DEFAULT_MIN_SW_SCORE "40"
#endif

#ifndef DEFAULT_MIN_SW_DIFF
#define DEFAULT_MIN_SW_DIFF "8"
#endif

#ifndef DEFAULT_SW_MATCH_SCORE
#define DEFAULT_SW_MATCH_SCORE "2"
#endif

#ifndef DEFAULT_SW_MISMATCH_PENALTY
#define DEFAULT_SW_MISMATCH_PENALTY "2"
#endif

#ifndef DEFAULT_SW_GAP_OPEN_PENALTY
#define DEFAULT_SW_GAP_OPEN_PENALTY "3"
#endif

#ifndef DEFAULT_SW_GAP_EXT_PENALTY
#define DEFAULT_SW_GAP_EXT_PENALTY "1"
#endif

using StringView = boost::string_ref;

namespace Utils {

std::vector<StringView> Tokenize(StringView line, char delimiter);

template <class T>
bool StringViewTo(const StringView& s, T& t) {
    std::istringstream ss(s.to_string());
    ss >> t;
    return ss.fail() ? false : true;
}

void Warning(const std::string& s);
void Error(const std::string& s);

}
#endif //SPLIT_POLYT_FROM_PBBAM_COMMON_HPP
