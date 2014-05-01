#include "creativity/belief/Profit.hpp"
#include <eris/algorithms.hpp>

using Eigen::RowVectorXd;

namespace creativity { namespace belief {

Profit::Profit(unsigned int D, Matrix<double, K, 1> beta_prior, double s_prior, Matrix<double, K, K> V_prior, double n_prior)
    : Linear<K>(std::move(beta_prior), std::move(s_prior), std::move(V_prior), std::move(n_prior)), D_{D}
{}

double Profit::predict(const double &q, const unsigned long &previousBooks, const unsigned long &marketBooks) {
    if (q < 0) throw std::domain_error("Profit::predict: q cannot be < 0");
    RowVectorXd X{K};
    X << 1, q, previousBooks == 0 ? 1 : 0, previousBooks, marketBooks;

    return X * beta_;
}

double Profit::argmax(
        const std::function<double(const double &)> q,
        const unsigned long &previousBooks, const unsigned long &marketBooks,
        const double &l_max) {
    const double Xb = beta_[0] + beta_[2] * (previousBooks == 0 ? 1 : 0) + beta_[3] * previousBooks + beta_[5] * marketBooks;
    auto profit = [&] (const double &l) -> double {
        return Xb + beta_[1] * q(l) - l;
    };

    return eris::single_peak_search(profit, 0, l_max);
}
        



}}
