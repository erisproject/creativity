#include <eris/algorithms.hpp>
#include <Eigen/Core>
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace eris;
using Eigen::RowVectorXd;
using Eigen::VectorXd;

const double alpha0 = -5, alpha1 = 3, alpha2 = 0.25;

double q(const double &l) {
    return alpha0 + alpha1*std::pow(l, alpha2);
}

int main() {
    double firstBook = 0, prevBooks = 3;
    double beta = 2;
    int d = 2;

    RowVectorXd Xt{3};
    VectorXd g{3};
    Xt << 1, firstBook, prevBooks;
    g << 10, -10, 2;

    auto profit = [&] (const double &l) -> double {
        double Nd = Xt * g + beta*q(l);
        return (Nd <= 0 ? 0 : std::pow(Nd, d)) - l;
    };

    std::cout << std::setprecision(10);
    std::cout << "profit(l) = [ (" << Xt << ") (" << g.transpose() << ")' + " << beta << " * ("
        << alpha0 << " + " << alpha1 << " l^" << alpha2 << ") ]^" << d << " - l\n";

    double max_l = single_peak_search(profit, 0, 1000);

    std::cout << "argmax = " << max_l << ", profit = " << profit(max_l) << "\n";
}
