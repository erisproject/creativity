#include "creativity/belief/ProfitStream.hpp"

using namespace eris;
using Eigen::RowVectorXd;

namespace creativity { namespace belief {

ProfitStream::ProfitStream(size_t K)
    : ProfitStream(VectorXd::Zero(K), 1.0, MatrixXd::Identity(K, K), 1e-6)
{
    beta_[K-1] = 1.0;
}

ProfitStream::ProfitStream(
        const Ref<const VectorKd> &beta,
        double s2,
        const Ref<const MatrixKd> &V,
        double n)
    : LinearBase(beta, s2, V, n)
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
