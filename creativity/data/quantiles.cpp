#include "creativity/data/quantiles.hpp"
#include <stdexcept>

namespace creativity { namespace data {

const std::regex
    quantile_field_regex("m(?:in|ax|edian)|q\\d*[1-9]"),
    ordinal_regex("(?:\\d*[02-9])?(?:1st|2nd|3rd)|\\d*1[123]th|\\d*[04-9]th");

std::string quantile_field(double quantile) {
    if (quantile == 0) return "min";
    else if (quantile == 0.5) return "median";
    else if (quantile == 1) return "max";
    else if (quantile > 0 and quantile < 1) return "q" + double_str(quantile, std::numeric_limits<double>::digits10).substr(2);
    
    throw std::invalid_argument("invalid argument passed to creativity::data::quantile_field: quantiles must be >= 0 and <= 1");
}

}}
