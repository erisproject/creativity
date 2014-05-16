#include "creativity/belief/Demand.hpp"
#include "creativity/Book.hpp"
#include "creativity/BookMarket.hpp"
#include "creativity/Reader.hpp"

using Eigen::RowVectorXd;

namespace creativity { namespace belief {

double Demand::predict(const double &P, const double &q, const unsigned long &S,
        const unsigned long &otherBooks, const unsigned long &marketBooks) const {
    if (P < 0) throw std::domain_error("Demand::predict: P cannot be < 0");
    if (q < 0) throw std::domain_error("Demand::predict: q cannot be < 0");
    RowVectorKd X;
    X << 1, std::pow(P, D_), std::pow(q, D_), S, otherBooks == 0 ? 1 : 0, otherBooks, marketBooks;

    return LinearBase::predict(X);
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


Demand::RowVectorKd Demand::bookRow(eris::SharedMember<Book> book, double quality) {
    RowVectorKd row{K()};

    auto t = book->simulation()->t() - 1;
    row << 1.0,
        (book->hasMarket() ? std::pow(book->market()->price(), D_) : 0.0),
        std::pow(quality, D_),
        book->lifeSales() - book->sales(t),
        book->age() - 1.0, // -1 because we're called in interopt, when t() has been incremented.
        book->author()->wrote().size() == 1 ? 1.0 : 0.0,
        book->author()->wrote().size() - 1,
        book->simulation()->countMarket<BookMarket>();



            MatrixXKd X(books.size(), K());
            size_t i = 0;
            for (const eris::SharedMember<Book> &book : books) {
                X(i, 0) = 1;
                X(i, 1) = (book->hasMarket() ? std::pow(book->market()->price(), D_) : 0);
                X(i, 2) = book->order();
                X(i, 3) = book->age();
                double price = book->hasMarket() ? book->market()->price() : 0.0;
                X(i, 4) = price;
                X(i, 5) = price*book->age();
                X(i, 6) = book->lifeSales();
                i++;
            }
            return X;
        }

}}
