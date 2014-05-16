#include "creativity/belief/Profit.hpp"
#include <eris/algorithms.hpp>

#include <eris/debug.hpp>

namespace creativity { namespace belief {

double Profit::predict(const double &q, const unsigned long &previousBooks, const unsigned long &marketBooks) const {
    if (q < 0) throw std::domain_error("Profit::predict: q cannot be < 0");
    RowVectorKd X{K()};
    X << 1, std::pow(q, D_), previousBooks == 0 ? 1 : 0, previousBooks, marketBooks;
    return LinearBase::predict(X);
}

double Profit::argmaxL(
        const std::function<double(const double &)> q,
        const unsigned long &previousBooks, const unsigned long &marketBooks,
        const double &l_max
) const {

    // Get the prediction without the quality term:
    const double Xb_base = predict(0, previousBooks, marketBooks);
    auto profit = [&] (const double &l) -> double {
        double quality = q(l);
        // Preserve the sign across the exponentiation
        return Xb_base + beta_[1] * std::copysign(std::pow(quality, D_), quality) - l;
    };

    return eris::single_peak_search(profit, 0, l_max);
}


}}
