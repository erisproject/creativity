#include "creativity/belief/Demand.hpp"
#include "creativity/Book.hpp"
#include "creativity/BookMarket.hpp"
#include "creativity/Reader.hpp"

using namespace Eigen;
using namespace eris;

namespace creativity { namespace belief {

unsigned int Demand::fixedModelSize() const { return parameters(); }

double Demand::predict(double P, double q, unsigned int S, unsigned int last_sales, unsigned int age, unsigned int otherBooks, unsigned int marketBooks) {
    if (P < 0) throw std::domain_error("Demand::predict: P cannot be < 0");
    if (q < 0) throw std::domain_error("Demand::predict: q cannot be < 0");
    RowVectorXd X(K());
    X << 1.0, P, q, q*q, S, (last_sales == 0 and age > 0) ? 1 : 0, age == 0 ? 1 : 0, age, otherBooks == 0 ? 1 : 0, otherBooks, marketBooks;

    return predict(X);
}

std::pair<double, double> Demand::argmaxP(
        double q, unsigned int s, unsigned int z, unsigned int n, unsigned int last_sales, unsigned int age,
        unsigned int otherBooks, unsigned int marketBooks, double c) {
    if (c < 0) throw std::domain_error("Demand::argmaxP: `c` must be non-negative");

    int demand_max = n - 1 - s - z;
    if (demand_max <= 0) return {0, 0}; // The whole world already has a copy!

    // Add up all the fixed parts by setting P = 0:
    const double Xg = predict(0, q, s, last_sales, age, otherBooks, marketBooks);

    if (Xg <= 0)
        // The prediction for a price of 0 is 0 or negative, and any positive price is just going to
        // make that more negative, so give up: there's no demand left.
        return {0, 0};

    const double b = mean_beta_[1]; // negative
    double pmax = Xg / (-2*b) + c/2;
    if (Xg + b*pmax > demand_max) {
        // Predicted demand exceeds the demand cut-off, so choose the price at the cutoff
        pmax = (Xg - demand_max) / (-b);
    }
    const double q_pred = Xg + b*pmax;

    if (pmax > 1000) ERIS_DBGVAR(mean_beta_.transpose());
    if (pmax <= c or q_pred <= 0) return {0, 0};

    return {pmax, q_pred};
}


RowVectorXd Demand::bookRow(SharedMember<Book> book, double quality, unsigned int lag_market_books) const {
    RowVectorXd row(K_);

    if (book->simulation()->runStage() != RunStage::inter_Optimize) throw std::logic_error("Demand::bookRow() must be called during inter-optimization stage");

    if (not book->hasMarket()) throw std::logic_error("Demand::bookRow() cannot be called for off-market book");

    auto t = book->simulation()->t() - 1;

    if (book->sales(t+1) > 0) throw std::logic_error("Demand::bookRow() called with sales in period t (this should be impossible: sales shouldn't have happened in interoptimization stages!)");

    row << 1.0, book->price(), quality, quality*quality,
        book->lifeSales() - book->sales(t), // sales before time t
        (book->age() > 1 and book->sales(t) == 0) ? 1.0 : 0.0,
        book->age() == 1 ? 1.0 : 0.0,
        book->age() - 1.0, // age (-1 because we're called after t() has been incremented)
        book->author()->wrote().size() == 1 ? 1.0 : 0.0, // onlyBook
        book->author()->wrote().size() - 1, // otherBooks
        lag_market_books; // marketBooks

    return row;
}

}}
