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
    X << 1, std::pow(P, D_), std::pow(q, D_), S, otherBooks == 0 ? 1 : 0, otherBooks, marketBooks;

    return X * beta_;
}

double Demand::argmaxP(const double &q, const unsigned long &S, const unsigned long &otherBooks, const unsigned long &marketBooks,
                const double &c) {
    if (c < 0) throw std::domain_error("Demand::argmaxP: `c` must be non-negative");
    if (q < 0) throw std::domain_error("Demand::argmaxP: `q` must be non-negative");

    const unsigned int onlyBook = otherBooks == 0 ? 1 : 0;
    // Add up all the fixed parts:
    const double Xg =
        beta_[0] +
        // Purposely missing beta_1 P^D
        beta_[2] * std::pow(q, D_) +
        beta_[3] * S +
        beta_[4] * onlyBook +
        beta_[5] * otherBooks +
        beta_[6] * marketBooks;

    if (Xg <= 0)
        // Our prediction is a 0 (or negative) quantity, so just give back c
        return c;

    if (c == 0) {
        // With c = 0 the solution is straighforward:
        double P_d = Xg / (-beta_[1] * (1 + D_));
        return std::pow(P_d, 1.0/D_);
    }

    // Otherwise optimize numerically
    return eris::single_peak_search([&c,&Xg,this] (const double &P) -> double {
            return (P - c) * (Xg + beta_[1] * std::pow(P, D_));
            }, c, argmaxP_MAX);
}

}}
