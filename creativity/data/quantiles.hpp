#pragma once
#include "creativity/data/util.hpp"
#include <string>
#include <limits>
#include <regex>

namespace creativity { namespace data {

/** Takes a quantile and returns a CSV header representing that quantile.  The special quantiles 0,
 * 0.5, and 1 become "min", "median", and "max", respectively.  All other values take the form qXXX
 * where XXX is formed by converting the quantile to a string (via creativity::data::double_str)
 * have the leading "0." and any trailing 0s removed.  Thus 0.052 becomes q052, 0.52 becomes q52,
 * and 0.999 becomes q999.
 *
 * \param quantile the quantile to convert to header form
 * \returns the header string version of the quantile
 * \sa quantile_field_regex
 * \throws std::invalid_argument if given a quantile not in [0,1]
 */
std::string quantile_field(double quantile);

/** A std::regex object that matches all values returned by quantile_field().
 * \sa quantile_field()
 */
extern const std::regex quantile_field_regex;

/** A std::regex object that matches ordinal values as produced by creativity-series. */
extern const std::regex ordinal_regex;


}}
