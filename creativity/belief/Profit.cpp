#include "creativity/belief/Profit.hpp"
#include "creativity/BookMarket.hpp"
#include <eris/algorithms.hpp>

#include <eris/debug.hpp>

namespace creativity { namespace belief {

using namespace Eigen;

unsigned int Profit::fixedModelSize() const { return parameters(); }

double Profit::predict(unsigned int draws, double q, unsigned long previousBooks, unsigned long marketBooks) {
    if (q < 0) throw std::domain_error("Profit::predict(): illegal negative quality value");
    return predict(profitRow(q, previousBooks, marketBooks), draws);
}

double Profit::argmaxL(
        unsigned int draws,
        const std::function<double(const double &)> q,
        unsigned long previousBooks, unsigned long marketBooks,
        double l_max
) {

    // Get the part of the prediction without the quality term:
    RowVectorXd X = profitRow(0.0, previousBooks, marketBooks);
    auto profit = [&q,&X,&draws,this] (const double &l) -> double {
        double quality = q(l);
        if (quality < 0.0) throw std::logic_error("Profit::argmaxL(): quality function returned invalid negative value");
        X[1] = quality;
        X[2] = quality*quality;
        return predict(X, draws) - l;
    };

    return eris::single_peak_search(profit, 0, l_max);
}

RowVectorXd Profit::profitRow(double quality, int, int lag_market_books) {
    RowVectorXd Xi(parameters());
    Xi << 1.0, quality, quality*quality, lag_market_books;
    return Xi;
}

}}
