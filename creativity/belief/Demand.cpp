#include <eris/Simulation.hpp>
#include "creativity/belief/Demand.hpp"
#include "creativity/Book.hpp"
#include "creativity/Reader.hpp"
#include <stdexcept>
#include <eris/algorithms.hpp>

using namespace Eigen;
using namespace eris;

namespace creativity { namespace belief {

unsigned int Demand::fixedModelSize() const { return parameters(); }

double Demand::predict(unsigned int draws, double P, double q, unsigned int S, unsigned int nosales, unsigned int age, unsigned int otherBooks, unsigned int marketBooks) {
    if (P < 0) throw std::domain_error("Demand::predict: P cannot be < 0");
    if (q < 0) throw std::domain_error("Demand::predict: q cannot be < 0");

    return predict(row(P, q, S, nosales, age, otherBooks, marketBooks), draws)[0];
}

std::pair<double, double> Demand::argmaxP(
        unsigned int draws,
        double q, unsigned int s, unsigned int z, unsigned int n, unsigned int nosales, unsigned int age,
        unsigned int otherBooks, unsigned int marketBooks, double c,
        double max_price) {
    if (c < 0) throw std::domain_error("Demand::argmaxP: `c` must be non-negative");

    double demand_max = n - 1 - s - z;
    if (demand_max <= 0) return {0, 0}; // The whole world already has a copy!

    // Get an X row within which we'll change the price to determine the profit-maximizing price
    RowVectorXd X = row(0, q, s, nosales, age, otherBooks, marketBooks);
    if (predict(X, draws)[0] <= 0)
        // The prediction at a price of 0 is 0 or negative, and any positive price is just going to
        // make that more negative, so give up: there's no demand left.
        return {0, 0};

    auto profit = [&X,&draws,&c,&demand_max,this] (const double &p) -> double {
        X[1] = p;
        double q = predict(X, draws)[0];
        if (q > demand_max) q = demand_max;
        return (p - c) * q;
    };

    double pstar = eris::single_peak_search(profit, c, max_price);
    X[1] = pstar;
    double qstar = predict(X, draws)[0];
    if (qstar > demand_max) qstar = demand_max;
    double pistar = (pstar - c)*qstar;

    if (pistar < 0) return {0, 0};
    return {pstar, qstar};
}

RowVectorXd Demand::row(double P, double q, unsigned int s, unsigned int nosales,
        unsigned int age, unsigned int otherBooks, unsigned int marketBooks) {
    if (false and otherBooks and marketBooks) {} // silence warning
    RowVectorXd row(parameters());

    row << 1.0, P, q, /* q*q,*/ s, nosales > age ? age : nosales; //, age == 0 ? 1 : 0, age, otherBooks, marketBooks;

    return row;
}

RowVectorXd Demand::bookRow(SharedMember<Book> book, double quality, unsigned int lag_market_books) {
    if (book->simulation()->runStage() != Simulation::RunStage::inter_Optimize) throw std::logic_error("Demand::bookRow() must be called during inter-optimization stage");

    if (not book->hasPrivateMarket()) throw std::logic_error("Demand::bookRow() cannot be called for off-market book");

    auto t = book->simulation()->t() - 1;

    if (book->sales(t+1) > 0) throw std::logic_error("Demand::bookRow() called with sales in period t (this should be impossible: sales shouldn't have happened in interoptimization stages!)");
    if (book->lastSale() > t) throw std::logic_error("Demand::bookRow() called with last sale in period > t (this should be impossible: sales shouldn't have happened in interoptimization stages!)");
    unsigned int nosales = t - book->lastSale();
    if (nosales == t) nosales = book->age();

    return row(book->price(), quality, book->lifeSales() - book->sales(t), nosales, book->age(), book->author()->wrote().size() - 1, lag_market_books);
}

}}
