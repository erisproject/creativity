#include "creativity/belief/ProfitStream.hpp"

using Eigen::RowVectorXd;
using namespace eris;

namespace creativity { namespace belief {

void ProfitStream::verifyParameters() const {
    for (unsigned int i = 0; i < K(); i++) {
        if (beta_[i] < 0) throw std::domain_error("ProfitStream beta_prior must be non-negative");
    }
}

double ProfitStream::predict(double profit_curr, unsigned int age) const {
    if (profit_curr <= 0 or age >= K()) return 0.0;
    return beta_[age] * profit_curr;
}

void ProfitStream::populate(SharedMember<Book> book, Ref<VectorXd> y, Ref<MatrixXKd> X, size_t n) const {
    double total_revenue = book->lifeRevenue();
    double cumul_revenue = 0;
    const unsigned long &t0 = book->created();
    const auto k = K();
    X.block(k*n, 0, k, k).setZero();
    y.segment(k*n, k).setZero();
    for (int i = 0; i < K(); i++) {
        cumul_revenue += book->revenue(t0 + i);
        double profit_remaining = total_revenue - cumul_revenue;
        // Check for numerical error
        if (profit_remaining < 0) profit_remaining = 0;

        y[k*n+i] = profit_remaining; // Profit remaining
        X(k*n+i, i) = cumul_revenue; // Profit so far
    }
}

ProfitStream ProfitStream::update(const Ref<const VectorXd> &y, const Ref<const MatrixXKd> &X) const {
    return ProfitStream{LinearBase::update(y, X)};
}

}}
