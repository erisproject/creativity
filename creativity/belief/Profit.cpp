#include "creativity/belief/Profit.hpp"
#include <eris/algorithms.hpp>
#include <stdexcept>

namespace creativity { namespace belief {

using namespace Eigen;

unsigned int Profit::fixedModelSize() const { return parameters(); }

double Profit::predict(unsigned int draws, double q, unsigned long previousBooks, unsigned long marketBooks) {
    if (q < 0) throw std::domain_error("Profit::predict(): illegal negative quality value");
    return predict(profitRow(q, previousBooks, marketBooks), draws)[0];
}

std::pair<double, double> Profit::argmaxL(
        unsigned int draws,
        const std::function<double(const double &)> q,
        unsigned long previousBooks,
        unsigned long marketBooks,
        double l_min,
        double l_max
) {

    if (l_min > l_max) throw std::logic_error("Profit::argmaxL called with l_min > l_max");

    // Get the part of the prediction without the quality term:
    RowVectorXd X = profitRow(0.0, previousBooks, marketBooks);
    auto profit = [&q,&X,&draws,this] (const double &l) -> double {
        double quality = q(l);
        if (quality < 0.0) throw std::logic_error("Profit::argmaxL(): quality function returned invalid negative value");
        X[1] = quality;
        X[2] = quality*quality;
        return predict(X, draws)[0] - l;
    };

    double ell = eris::single_peak_search(profit, l_min, l_max);
    double net_profit = predict(draws, q(ell), previousBooks, marketBooks) - ell;
    return std::pair<double, double>(ell, net_profit);
}

RowVectorXd Profit::profitRow(double quality, int /*previous_books*/, int lag_market_books) {
    RowVectorXd Xi(parameters());
    if (quality < 0) quality = 0;
    Xi << 1.0, quality, quality*quality, lag_market_books;
    return Xi;
}

}}
