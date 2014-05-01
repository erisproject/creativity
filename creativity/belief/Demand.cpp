#include "creativity/belief/Demand.hpp"

using Eigen::RowVectorXd;

namespace creativity { namespace belief {

Demand::Demand(unsigned int D, Matrix<double, K, 1> beta_prior, double s_prior, Matrix<double, K, K> V_prior, double n_prior)
    : Linear<K>(std::move(beta_prior), std::move(s_prior), std::move(V_prior), std::move(n_prior)), D_{D}
{}

double Demand::predict(const double &P, const double &q, const unsigned long &S,
        const unsigned long &otherBooks, const unsigned long &marketBooks) {
    if (P < 0) throw std::domain_error("Demand::predict: P cannot be < 0");
    if (q < 0) throw std::domain_error("Demand::predict: q cannot be < 0");
    RowVectorXd X{K};
    X << 1, P, q, S, otherBooks == 0 ? 1 : 0, otherBooks, marketBooks;

    return X * beta_;
}

}}
