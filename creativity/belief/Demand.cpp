#include "creativity/belief/Demand.hpp"
#include "creativity/Book.hpp"
#include "creativity/BookMarket.hpp"
#include "creativity/Reader.hpp"

using Eigen::RowVectorXd;

namespace creativity { namespace belief {

double Demand::predict(double P, double q, unsigned long S, unsigned long otherBooks, unsigned long marketBooks) const {
    if (P < 0) throw std::domain_error("Demand::predict: P cannot be < 0");
    RowVectorKd X;
    X << 1, std::pow(P, D_), std::copysign(std::pow(q, D_), q), S, 0, otherBooks == 0 ? 1 : 0, otherBooks, marketBooks;

    return LinearBase::predict(X);
}

std::pair<double, double> Demand::argmaxP(double q, unsigned long S, unsigned long otherBooks, unsigned long marketBooks, double c) const {
    if (c < 0) throw std::domain_error("Demand::argmaxP: `c` must be non-negative");

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

    ERIS_DBGVAR(beta_[1]);

    // Otherwise optimize numerically
    double p_opt = eris::single_peak_search([&c,&Xg,this] (double P) -> double {
            return (P - c) * (Xg + beta_[1] * std::pow(P, D_));
            }, c, argmaxP_MAX);
    return {p_opt, Xg + beta_[1] * std::pow(p_opt, D_)};
}


Demand::RowVectorKd Demand::bookRow(eris::SharedMember<Book> book, double quality) const {
    RowVectorKd row{K()};

    ERIS_DBGVAR(K());

    auto t = book->simulation()->t() - 1;
    row << 1.0,
        (book->hasMarket() ? std::pow(book->market()->price(), D_) : 0.0), // Price^D
        std::pow(quality, D_), // quality^D
        book->lifeSales() - book->sales(t), // sales before time t
        book->age() - 1.0, // age (-1 because we're called agter t() has been incremented)
        book->author()->wrote().size() == 1 ? 1.0 : 0.0, // onlyBook
        book->author()->wrote().size() - 1, // previousBooks
        book->simulation()->countMarkets<BookMarket>(); // marketBooks

    return row;
}

Demand Demand::update(const Ref<const VectorXd> &y, const Ref<const MatrixXKd> &X) const {
    ERIS_DBG("update");
    return Demand{D_, LinearBase::update(y, X)};
}

}}
