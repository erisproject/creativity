#include "creativity/belief/Profit.hpp"
#include "creativity/BookMarket.hpp"
#include <eris/algorithms.hpp>

#include <eris/debug.hpp>

namespace creativity { namespace belief {

using namespace Eigen;

unsigned int Profit::fixedModelSize() const { return parameters(); }

double Profit::predict(double q, unsigned long previousBooks, unsigned long marketBooks) {
    if (q < 0) throw std::domain_error("Profit::predict(): illegal negative quality value");
    RowVectorXd X(K_);
    X << 1.0, q, q*q, previousBooks == 0 ? 1.0 : 0.0, previousBooks, marketBooks;

    return predict(X);
}

double Profit::argmaxL(
        const std::function<double(const double &)> q,
        unsigned long previousBooks, unsigned long marketBooks,
        double l_max
) {

    // Get the part of the prediction without the quality term:
    const double Xb_base = predict(0.0, previousBooks, marketBooks);
    auto profit = [&q,&Xb_base,this] (const double &l) -> double {
        double quality = q(l);
        if (quality < 0.0) throw std::logic_error("Profit::argmaxL(): quality function returned invalid negative value");
        return Xb_base + beta_[1] * quality + beta_[2] * quality * quality - l;
    };

    return eris::single_peak_search(profit, 0, l_max);
}

RowVectorXd Profit::profitRow(eris::SharedMember<Book> book, double quality) const {
    double prev_books = book->order();
    double first_book = prev_books == 0 ? 1.0 : 0.0;
    double market_books = book->simulation()->countMarkets<BookMarket>();

    RowVectorXd Xi(K());
    Xi << 1.0, quality, quality*quality, first_book, prev_books, market_books;
    return Xi;
}

}}
