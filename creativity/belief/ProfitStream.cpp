#include "creativity/belief/ProfitStream.hpp"

using namespace eris;
using Eigen::RowVectorXd;

namespace creativity { namespace belief {

using namespace Eigen;

double ProfitStream::predict(SharedMember<Book> book) {
    RowVectorXd X(K());
    for (size_t i = 0; i < K(); i++) {
        X[i] = book->revenue(book->created() + i);
    }

    return Linear::predict(X);
}

ProfitStream ProfitStream::newDerived(Linear &&model) const {
    return ProfitStream(std::move(model));
}



} }
