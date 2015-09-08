#include "creativity/cmdargs/strings.hpp"
#include <regex>
#include <string>
#include <sstream>

namespace creativity { namespace cmdargs {

template <>
std::string output_string(double v) {
    return std::regex_replace(
            std::regex_replace(std::to_string(v),
                std::regex("(\\.\\d*?)0+$"),
                "$1"),
            std::regex("\\.$"),
            "");
}

}}
