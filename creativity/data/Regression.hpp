#pragma once
#include <stdexcept>

namespace creativity { namespace data {

/** Exception class for rank errors, for example when trying to do OLS with n < k, or an X matrix
 * without full column rank.
 */
class RankError : public std::runtime_error {
    public:
        RankError(const std::string &message) : std::runtime_error(message) {}
};

}}
