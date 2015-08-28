#pragma once
#include <string>
#include <type_traits>

namespace creativity { namespace cmdargs {

/** Returns an argument name for T.  For numeric types, this is one of 𝑹 (for floating point
 * values), 𝑵 (for unsigned integer types), or 𝒁 (for signed integer types).  For other
 * types, this is simply "arg".
 */
template <typename T>
std::string type_string() {
    return std::is_floating_point<T>::value ? u8"ℝ" :
        std::is_integral<T>::value ? std::is_unsigned<T>::value ? u8"ℕ" : u8"ℤ" :
        u8"arg";
}

/** Returns a value converted to a string; in most cases this just passes the value to
 * std::to_string for conversion.
 *
 * A specialization of this for doubles converts the double to a string via
 * std::to_string(), then trims off trailing 0's and, possibly, a trailing `.`.
 */
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
std::string output_string(T v) {
    return std::to_string(v);
}

/** Specialization of output_string for double types that trims trailing 0's and a trailing decimal
 * point from the string representation before returning it.
 */
template <> std::string output_string(double v); // Specialization in .cpp

}}
