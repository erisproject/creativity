#include "creativity/belief/Demand.hpp"
#include "creativity/Book.hpp"
#include "creativity/BookMarket.hpp"
#include "creativity/Reader.hpp"

using Eigen::RowVectorXd;

namespace creativity { namespace belief {

using namespace Eigen;

unsigned int Demand::fixedModelSize() const { return 8; }

double Demand::predict(double P, double q, unsigned long S, unsigned long otherBooks, unsigned long marketBooks) const {
    ERIS_DBGVAR(P);
    if (P < 0) throw std::domain_error("Demand::predict: P cannot be < 0");
    RowVectorXd X(K());
    X << 1, std::pow(P, D_), std::copysign(std::pow(q, D_), q), S, 0, otherBooks == 0 ? 1 : 0, otherBooks, marketBooks;

    ERIS_DBG("hi");
    return Linear::predict(X);
}

std::pair<double, double> Demand::argmaxP(double q, unsigned long S, unsigned long otherBooks, unsigned long marketBooks, double c) const {
    ERIS_DBG("hi");
    if (c < 0) throw std::domain_error("Demand::argmaxP: `c` must be non-negative");

    ERIS_DBG("hi");
    // Add up all the fixed parts by setting P = 0:
    const double Xg = predict(0, q, S, otherBooks, marketBooks);

    ERIS_DBGVAR(Xg);
    if (Xg <= 0)
        // Our prediction with a price of 0 is 0 (or negative) quantity, and any positive price is
        // just going to make that more negative, so just give back c
        return {c, Xg};

    ERIS_DBGVAR(c);
    if (c == 0) {
        // With c = 0 the solution is analytical:
        double P_d = Xg / (-beta_[1] * (1 + D_));
        return {std::pow(P_d, 1.0/D_), Xg + beta_[1] * P_d};
    }

    ERIS_DBG("hi");
    // Otherwise optimize numerically
    double p_opt = eris::single_peak_search([&c,&Xg,this] (double P) -> double {
            return (P - c) * (Xg + beta_[1] * std::pow(P, D_));
            }, c, argmaxP_MAX);
    ERIS_DBGVAR(p_opt);
    return {p_opt, Xg + beta_[1] * std::pow(p_opt, D_)};
}


RowVectorXd Demand::bookRow(eris::SharedMember<Book> book, double quality) const {
    RowVectorXd row(K());

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

Demand Demand::update(const Ref<const VectorXd> &y, const Ref<const MatrixXd> &X) const {
    return Demand{D_, Linear::update(y, X)};
}

}}
