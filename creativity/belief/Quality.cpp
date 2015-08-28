#include "creativity/belief/Quality.hpp"
#include "creativity/Book.hpp"

using Eigen::RowVectorXd;

namespace creativity { namespace belief {

using namespace Eigen;

unsigned int Quality::fixedModelSize() const { return parameters(); }

double Quality::predict(const Book &book, unsigned int draws) {
    RowVectorXd X(K());
    X << book.order(), book.age(), book.lifeSales() + book.lifePirated();
    return predict(X, draws);
}

}}
