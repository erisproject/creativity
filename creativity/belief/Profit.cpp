#include "creativity/belief/Profit.hpp"
#include <eris/algorithms.hpp>

namespace creativity { namespace belief {

Profit::Profit(
        const unsigned int &D,
        const VectorKd &beta_prior,
        const double &s_prior,
        const MatrixKd &V_prior,
        const double &n_prior
        )
    : Linear<KK>(beta_prior, s_prior, V_prior, n_prior), D_{D}
{}

double Profit::predict(const double &q, const unsigned long &previousBooks, const unsigned long &marketBooks) const {
    if (q < 0) throw std::domain_error("Profit::predict: q cannot be < 0");
    RowVectorKd X;
    X << 1, std::pow(q, D_), previousBooks == 0 ? 1 : 0, previousBooks, marketBooks;
    return Linear<KK>::predict(X);
}

double Profit::argmaxL(
        const std::function<double(const double &)> q,
        const unsigned long &previousBooks, const unsigned long &marketBooks,
        const double &l_max) const {

    // Get the prediction without the quality term:
    const double Xb_base = predict(0, previousBooks, marketBooks);
    auto profit = [&] (const double &l) -> double {
        return Xb_base + beta_[1] * std::pow(q(l), D_) - l;
    };

    return eris::single_peak_search(profit, 0, l_max);
}


}}
