#include "creativity/belief/Quality.hpp"
#include "creativity/Book.hpp"
#include "creativity/BookMarket.hpp"

using Eigen::RowVectorXd;

namespace creativity { namespace belief {

using namespace Eigen;

unsigned int Quality::fixedModelSize() const { return parameters(); }

double Quality::predict(const Book &book) {
    double price = 0;
    if (book.hasMarket())
        price = book.market()->price();
    RowVectorXd X(K());
    X << 1, book.order() == 0, book.order(), book.age(), price, price*book.age(), book.lifeSales();
    return Linear::predict(X);
}


Quality Quality::newDerived(Linear &&model) const {
    return Quality(std::move(model));
}


}}
