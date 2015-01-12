#include "creativity/belief/Demand.hpp"
#include "creativity/Book.hpp"
#include "creativity/BookMarket.hpp"
#include "creativity/Reader.hpp"

using namespace Eigen;
using namespace eris;

namespace creativity { namespace belief {

unsigned int Demand::fixedModelSize() const { return parameters(); }

double Demand::predict(double P, double q, unsigned long S, unsigned long age, unsigned long otherBooks, unsigned long marketBooks) {
    if (P < 0) throw std::domain_error("Demand::predict: P cannot be < 0");
    if (q < 0) throw std::domain_error("Demand::predict: q cannot be < 0");
    RowVectorXd X(K());
    X << 1.0, P, q, q*q, S, age, otherBooks == 0 ? 1 : 0, otherBooks, marketBooks;

    return predict(X);
}

std::pair<double, double> Demand::argmaxP(double q, unsigned long S, unsigned long age, unsigned long otherBooks, unsigned long marketBooks, double c) {
    if (c < 0) throw std::domain_error("Demand::argmaxP: `c` must be non-negative");

    // Add up all the fixed parts by setting P = 0:
    const double Xg = predict(0, q, S, age, otherBooks, marketBooks);

    if (Xg < 0)
        // The prediction for a price of 0 is negative quantity, and any positive price is just
        // going to make that more negative, so give up: there's no demand left.
        return {0, 0};

    const double b = mean_beta_[1]; // negative
    const double pmax = Xg / (-2*b) + c/2;
    const double profit = (Xg + b*pmax)*(pmax - c);
    if (pmax < c or profit < 0) return {0, 0};

    return {pmax, profit};
}


RowVectorXd Demand::bookRow(SharedMember<Book> book, double quality) const {
    RowVectorXd row(K_);

    if (book->simulation()->runStage() != RunStage::inter_Optimize) throw std::logic_error("Demand::bookRow() must be called during inter-optimization stage");

    if (not book->hasMarket()) throw std::logic_error("Demand::bookRow() cannot be called for off-market book");

    auto t = book->simulation()->t() - 1;

    if (book->sales(t+1) > 0) throw std::logic_error("Demand::bookRow() called with sales in period t (this should be impossible: sales shouldn't have happened in interoptimization stages!)");

    row << 1.0, book->price(), quality, quality*quality,
        book->lifeSales() - book->sales(t), // sales before time t
        book->age() - 1.0, // age (-1 because we're called after t() has been incremented)
        book->author()->wrote().size() == 1 ? 1.0 : 0.0, // onlyBook
        book->author()->wrote().size() - 1, // otherBooks
        book->simulation()->countMarkets<BookMarket>(); // marketBooks

    return row;
}

Demand Demand::newDerived(Linear &&base) const {
    return Demand(D_, std::move(base));
}

}}
