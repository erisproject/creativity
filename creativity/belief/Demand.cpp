#include "creativity/belief/Demand.hpp"

using Eigen::RowVectorXd;

namespace creativity { namespace belief {

Demand::Demand(const unsigned int &D, const VectorKd &beta_prior, const double &s_prior, const MatrixKd &V_prior, const double &n_prior)
    : Linear<KK>(beta_prior, s_prior, V_prior, n_prior), D_{D}
{}

double Demand::predict(const double &P, const double &q, const unsigned long &S,
        const unsigned long &otherBooks, const unsigned long &marketBooks) const {
    if (P < 0) throw std::domain_error("Demand::predict: P cannot be < 0");
    if (q < 0) throw std::domain_error("Demand::predict: q cannot be < 0");
    Matrix<double, 1, KK> X;
    X << 1, std::pow(P, D_), std::pow(q, D_), S, otherBooks == 0 ? 1 : 0, otherBooks, marketBooks;

    return Linear<KK>::predict(X);
}

std::pair<double, double> Demand::argmaxP(const double &q, const unsigned long &S, const unsigned long &otherBooks, const unsigned long &marketBooks,
                const double &c) const {
    if (c < 0) throw std::domain_error("Demand::argmaxP: `c` must be non-negative");
    if (q < 0) throw std::domain_error("Demand::argmaxP: `q` must be non-negative");

    // Add up all the fixed parts by setting P = 0:
    const double Xg = predict(0, q, S, otherBooks, marketBooks);

    if (Xg <= 0)
        // Our prediction with a price of 0 is 0 (or negative) quantity, and any positive price is
        // just going to make that more negative, so just give back c
        return {c, Xg};

    if (c == 0) {
        // With c = 0 the solution is analytical:
        double P_d = Xg / (-beta_[1] * (1 + D_));
        return {std::pow(P_d, 1.0/D_), Xg + beta_[1] * P_d};
    }

    // Otherwise optimize numerically
    double p_opt = eris::single_peak_search([&c,&Xg,this] (const double &P) -> double {
            return (P - c) * (Xg + beta_[1] * std::pow(P, D_));
            }, c, argmaxP_MAX);
    return {p_opt, Xg + beta_[1] * std::pow(p_opt, D_)};
}

}}
