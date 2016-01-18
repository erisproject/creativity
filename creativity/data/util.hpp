#pragma once
#include <limits>
#include <string>
#include <sstream>
#include <vector>
#include <Eigen/Core>

/** \file
 * Miscellaneous utility functions for dealing with data.
 */

namespace creativity { namespace data {

/** Converts a double value to a string, use as much precision is required (but no more).  When
 * outputting values, we want maximum precision, but don't necessarily need the full max_digits10
 * digits to get there.  For example, 0.1 with max_digits10 significant digits becomes
 * 0.10000000000000001 but 0.1 also converts to the (numerically) identical value.
 * 0.10000000000000002, on the other hand, is a numerically distinct value and thus needs every
 * decimal digit.
 *
 * This function first tries to convert the value at the requested precision, then at the requested
 * precision less one, then less two; if the subsequent values are numerically identical to the
 * given double value when reconverted to a double, the shortest is returned.
 *
 * \param d the double to convert to a string
 * \param precision the requested precision, defaulting to std::numeric_limits<double>::max_digits10
 */
std::string double_str(double d, unsigned precision = std::numeric_limits<double>::max_digits10);

/** Returns a sample quantile from the vector of values.  Linear interpolation is used for quantiles
 * that lie between elements.  The interpolation is such that the quantile for 0 returns the
 * smallest value and the quantile for 1 returns the largest value.
 *
 * \param vals the vector of pre-sorted double values
 * \param prob the probability for which the quantile should be calculated; must be between 0 and 1.
 * \returns NaN if the vector is empty, otherwise returns the requested quantile
 */
double quantile(const std::vector<double> &vals, double prob);

/** Same as above, but operates on an Eigen vector-like object. */
double quantile(const Eigen::Ref<const Eigen::VectorXd> &vals, double prob);

/** Calculates and returns the sample variance of the given vector of values.  If the mean of the
 * vector is already known it can be passed in as mean, otherwise the default, NaN, calculates the
 * mean from the vector.
 */
double variance(const std::vector<double> &vals, double mean = std::numeric_limits<double>::quiet_NaN());

/** enum of the different strategies for dealing with prefix/suffix overlap in common_ends(). */
enum class OverlapReducer {
    None, ///< Don't reduce the prefix/suffix: the returned prefix and suffix may overlap
    Prefix, ///< Reduce the prefix by the minimum amount needed to eliminate the overlap
    Suffix, ///< Reduce the suffix by the minimum amount needed to eliminate the overlap
    Equalize, ///< Apply reductions to the prefix/suffix values that make the values closer to equal
    Disequalize ///< Apply reductions to whichever of the prefix/suffix is smaller
};

/** Determines the number of equal (by ==) front and/or back elements in a range of random-access
 * iterable objects.  Returns a pair where `.first` is the number of common elements from the
 * beginning of each given container, and `.second` is the number of common elements from the end of
 * each given container.
 *
 * For example:
 *
 *     std::list<std::string> strings({"abcdxyz", "abcxyz", "abcsomethingyz"});
 *     common_ends(strings); // Returns std::pair<size_t,size_t>(3,2)
 *
 * In the case of overlaps, for example with the strings "aabbc" and "aabbbaabbc" (which each have a
 * common prefix of length 4 ("aabb") and common suffix of length 5 ("aabbc")), the optional
 * overlap_reduce parameter controls how to reduce the overlap.
 *
 * \param it an iterator to the first element to consider
 * \param last the past-the-end iterator at which to stop
 * \param prefix whether to look for a common prefix (defaults to true).  If false, the returned
 * prefix length will be 0.
 * \param suffix whether to look for a common suffix (defaults to true).  If false, the returned
 * suffix length will be 0.
 * \param reducer an OverlapReducer value controlling what happens when the prefix and suffix
 * overlaps for one or more elements.  The default, OverlapReducer::Disequalize, reduces reduces
 * whichever is smaller (leaving the larger match unadjuated).
 * by approximately equal values to eliminate the overlap.
 * \returns a pair of size_t values where `.first` is the common prefix length and `.second` is
 * the common suffix length.
 */
template <typename InputIt>
typename std::enable_if<
    // InputIt must be an input iterator:
    std::is_base_of<std::input_iterator_tag, typename std::iterator_traits<InputIt>::iterator_category>::value
    and
    // The value_type of InputIt must itself be bidirectional-iterable
    std::is_base_of<
        std::bidirectional_iterator_tag,
        typename std::iterator_traits<typename InputIt::value_type::iterator>::iterator_category
    >::value,
    std::pair<size_t, size_t>
>::type
common_ends(InputIt it, InputIt last,
        bool prefix = true, bool suffix = true,
        OverlapReducer reducer = OverlapReducer::Disequalize) {
    if (it == last) return {0,0};
    const auto &ref_0 = *it;
    std::pair<size_t, size_t> common(0, 0);
    // Use the first element as our reference element, using its length as the initial
    // prefix/suffix values.
    auto shortest = std::distance(it->begin(), it->end());
    if (prefix) common.first = shortest;
    if (suffix) common.second = shortest;

    // Now, starting from the second element, check how much of the prefix/suffix matches
    for (it++; it != last and (common.first > 0 or common.second > 0); it++) {
        if (common.first > 0) {
            // Iterate from the beginning, and keep iterating as long as each element is equal along
            // the way
            auto ref_it = ref_0.begin();
            auto new_it = it->begin();
            // NB: don't need to check ref_it != ref_0.end() here (because the i < common.first will
            // ensure we can't go past the end).
            for (size_t i = 0; i < common.first and new_it != it->end(); i++, ref_it++, new_it++) {
                if (not (*new_it == *ref_it)) {
                    // We found two unequal elements at [i], so reduce the common prefix length to i
                    common.first = i;
                    break;
                }
            }
        }
        if (common.second > 0) {
            // Iterate backwards from the end
            auto ref_it = ref_0.rbegin();
            auto new_it = it->rbegin();
            for (size_t i = 0; i < common.second and new_it != it->rend(); i++, ref_it++, new_it++) {
                // We found two unequal elements at i from the end, so reduce the common suffix
                // length to i
                if (not (*new_it == *ref_it)) {
                    common.second = i;
                    break;
                }
            }
        }

        if (common.first > 0 and common.second > 0) {
            shortest = std::min(shortest, std::distance(it->begin(), it->end()));
        }
    }

    if (common.first + common.second > (size_t) shortest) {
        switch (reducer) {
            case OverlapReducer::None:
                break;
            case OverlapReducer::Prefix:
                common.first = shortest - common.second;
                break;
            case OverlapReducer::Suffix:
                common.second = shortest - common.first;
                break;
            case OverlapReducer::Equalize:
                while (common.first + common.second > (size_t) shortest) {
                    if (common.second >= common.first) common.second--;
                    else common.first--;
                }
                break;
            case OverlapReducer::Disequalize:
                if (common.second <= common.first) common.second = shortest - common.first;
                else common.first = shortest - common.second;
                break;
        }
    }

    return common;
}

}}
