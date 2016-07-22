#ifndef SPLIT_POLYT_FROM_PBBAM_COMMON_HPP
#define SPLIT_POLYT_FROM_PBBAM_COMMON_HPP

#include <vector>
#include <sstream>
#include <boost/utility/string_ref.hpp>

using StringView = boost::string_ref;

namespace Utils {

std::vector<StringView> Tokenize(StringView line, char delimiter);

template <class T>
bool StringViewTo(const StringView& s, T& t) {
    std::istringstream ss(s.to_string());
    ss >> t;
    return ss.fail() ? false : true;
}

}
#endif //SPLIT_POLYT_FROM_PBBAM_COMMON_HPP
