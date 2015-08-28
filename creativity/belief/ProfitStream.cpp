#include "creativity/belief/ProfitStream.hpp"
#include "creativity/Book.hpp"
#include <Eigen/Core>
#include <cstddef>

using namespace eris;
using Eigen::RowVectorXd;

namespace creativity { namespace belief {

using namespace Eigen;

double ProfitStream::predict(SharedMember<Book> book, unsigned int draws) {
    RowVectorXd X(K());
    for (size_t i = 0; i < K(); i++) {
        X[i] = book->revenue(book->created() + i);
    }

    return predict(X, draws);
}

} }
