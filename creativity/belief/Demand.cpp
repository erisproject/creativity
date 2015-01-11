#include "creativity/belief/Demand.hpp"
#include "creativity/Book.hpp"
#include "creativity/BookMarket.hpp"
#include "creativity/Reader.hpp"

using Eigen::RowVectorXd;

namespace creativity { namespace belief {

using namespace Eigen;

unsigned int Demand::fixedModelSize() const { return parameters(); }

double Demand::predict(double P, double q, unsigned long S, unsigned long otherBooks, unsigned long marketBooks) {
    if (P < 0) throw std::domain_error("Demand::predict: P cannot be < 0");
    if (q < 0) throw std::domain_error("Demand::predict: q cannot be < 0");
    RowVectorXd X(K());
    X << 1.0, P, P*P, q, q*q, S, 0, otherBooks == 0 ? 1 : 0, otherBooks, marketBooks;

    return predict(X);
}

std::pair<double, double> Demand::argmaxP(double q, unsigned long S, unsigned long otherBooks, unsigned long marketBooks, double c) {
    if (c < 0) throw std::domain_error("Demand::argmaxP: `c` must be non-negative");

    // Add up all the fixed parts by setting P = 0:
    const double Xg = predict(0, q, S, otherBooks, marketBooks);

    ERIS_DBGVAR(Xg);

    if (Xg <= 0)
        // Our prediction with a price of 0 is 0 (or negative) quantity, and any positive price is
        // just going to make that more negative, so just give back c
        return {c, Xg};

    // Otherwise optimize numerically
    double p_opt = eris::single_peak_search([&c,&Xg,this] (double P) -> double {
            return (P - c) * (Xg + beta_[1] * P + beta_[2] * P * P);
            }, c, argmaxP_MAX);
    return {p_opt, Xg + beta_[1] * p_opt + beta_[2] * p_opt * p_opt};
}


RowVectorXd Demand::bookRow(eris::SharedMember<Book> book, double quality) const {
    RowVectorXd row(K_);

    auto t = book->simulation()->t() - 1;
    double p = book->hasMarket() ? book->market()->price() : 0.0;
    row << 1.0, p, p*p, quality, quality*quality,
        book->lifeSales() - book->sales(t), // sales before time t
        book->age() - 1.0, // age (-1 because we're called agter t() has been incremented)
        book->author()->wrote().size() == 1 ? 1.0 : 0.0, // onlyBook
        book->author()->wrote().size() - 1, // previousBooks
        book->simulation()->countMarkets<BookMarket>(); // marketBooks

    return row;
}

Demand Demand::newDerived(Linear &&base) const {
    return Demand(D_, std::move(base));
}

}}
