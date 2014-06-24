#include "creativity/belief/ProfitStream.hpp"

using namespace eris;
using Eigen::RowVectorXd;

namespace creativity { namespace belief {

ProfitStream::ProfitStream(size_t K)
    : LinearBase(VectorXd::Zero(K), 1e20, 1e10 * MatrixXd::Identity(K, K), 1e-6)
{}

double ProfitStream::predict(SharedMember<Book> book) const {
    RowVectorXd X(K());
    for (size_t i = 0; i < K(); i++) {
        X[i] = book->revenue(book->created() + i);
    }

    return LinearBase::predict(X);
}

ProfitStream ProfitStream::update(const Ref<const VectorXd> &y, const Ref<const MatrixXKd> &X) const {
    return ProfitStream{LinearBase::update(y, X)};
}


} }
