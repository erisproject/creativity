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

}}
