#include "creativity/belief/Quality.hpp"
#include "creativity/Book.hpp"
#include "creativity/BookMarket.hpp"

using Eigen::RowVectorXd;

namespace creativity { namespace belief {

Quality::Quality(
                const VectorKd &beta_prior,
                const double &s_prior,
                const MatrixKd &V_prior,
                const double &n_prior
              )
    : Linear<KK>(beta_prior, s_prior, V_prior, n_prior)
{}

double Quality::predict(const Book &book) const {
    RowVectorXd X{7};
    double price = 0;
    if (book.hasMarket())
        price = book.market()->price();
    X << 1, book.order() == 0, book.order(), book.age(), price, price*book.age(), book.sales();
    return Linear<KK>::predict(X);
}


}}
