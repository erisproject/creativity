#include "creativity/belief/ProfitStream.hpp"

using namespace eris;
using Eigen::RowVectorXd;

namespace creativity { namespace belief {

using namespace Eigen;

ProfitStream::ProfitStream(unsigned int K)
    : ProfitStream(VectorXd::Zero(K), 1.0, MatrixXd::Identity(K, K), 1e-6)
{
    if (K == 0) throw std::domain_error("Unable to create ProfitStream belief with K=0 parameters");
    beta_[K-1] = 1.0;
}

ProfitStream::ProfitStream(
        const Ref<const VectorXd> &beta,
        double s2,
        const Ref<const MatrixXd> &V,
        double n)
    : Linear(beta, s2, V, n)
{}

double ProfitStream::predict(SharedMember<Book> book) {
    RowVectorXd X(K());
    for (size_t i = 0; i < K(); i++) {
        X[i] = book->revenue(book->created() + i);
    }

    return Linear::predict(X);
}

ProfitStream ProfitStream::update(const Ref<const VectorXd> &y, const Ref<const MatrixXd> &X) const {
    return ProfitStream{Linear::update(y, X)};
}


} }
