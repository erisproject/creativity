#include "creativity/belief/Profit.hpp"
#include "creativity/BookMarket.hpp"
#include <eris/algorithms.hpp>

#include <eris/debug.hpp>

namespace creativity { namespace belief {

using namespace Eigen;

unsigned int Profit::fixedModelSize() const { return 5; }

double Profit::predict(double q, unsigned long previousBooks, unsigned long marketBooks) const {
    RowVectorXd X(K());
    X << 1, std::copysign(std::pow(q, D_), q), previousBooks == 0 ? 1 : 0, previousBooks, marketBooks;
    return Linear::predict(X);
}

double Profit::argmaxL(
        const std::function<double(const double &)> q,
        unsigned long previousBooks, unsigned long marketBooks,
        double l_max
) const {

    // Get the prediction without the quality term:
    const double Xb_base = predict(0, previousBooks, marketBooks);
    auto profit = [&q,&Xb_base,this] (const double &l) -> double {
        double quality = q(l);
        // Preserve the sign across the exponentiation
        return Xb_base + beta_[1] * std::copysign(std::pow(quality, D_), quality) - l;
    };

    return eris::single_peak_search(profit, 0, l_max);
}

RowVectorXd Profit::profitRow(eris::SharedMember<Book> book, double quality) const {
    double prev_books = book->order();
    double first_book = prev_books == 0 ? 1 : 0;
    double market_books = book->simulation()->countMarkets<BookMarket>();

    RowVectorXd Xi(K());
    Xi << 1, std::pow(quality, D_), first_book, prev_books, market_books;
    return Xi;
}

Profit Profit::update(const Ref<const VectorXd> &y, const Ref<const MatrixXd> &X) const {
    return Profit{D_, Linear::update(y, X)};
}


}}
