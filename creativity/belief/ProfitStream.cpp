#include "creativity/belief/ProfitStream.hpp"

using Eigen::RowVectorXd;

namespace creativity { namespace belief {

void ProfitStream::verifyParameters() const {
    for (unsigned int i = 0; i < K; i++) {
        if (beta_[i] < 0) throw std::domain_error("ProfitStream beta_prior must be non-negative");
    }
}

double ProfitStream::predict(const double &profit_curr, const unsigned int &age) const {
    if (profit_curr <= 0 or age >= K) return 0.0;
    return beta_[age] * profit_curr;
}

}}
