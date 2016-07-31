#ifndef SPLIT_POLYT_FROM_PBBAM_COMMON_HPP
#define SPLIT_POLYT_FROM_PBBAM_COMMON_HPP

#include <vector>
#include <sstream>
#include <boost/utility/string_ref.hpp>
#include "kernel_color.hpp"


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
