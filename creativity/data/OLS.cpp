#include "creativity/data/OLS.hpp"
#include "creativity/data/Variable.hpp"
#include "creativity/data/Regression.hpp"
#include <Eigen/SVD>
#include <algorithm>
#include <limits>
#include <string>
#include <ostream>

using namespace Eigen;
namespace creativity { namespace data {

OLS::OLS(const Equation &model) : model_(model) {}
OLS::OLS(Equation &&model) : model_(std::move(model)) {}

const Equation& OLS::model() const { return model_; }

void OLS::gather() {
    if (gathered_) return;

    y_ = model_.depVar().values();
    X_ = MatrixXd(y_.size(), model_.numVars());

    unsigned int k = 0;
    for (const Variable &var : model_) {
        var.populate(X_.col(k++));
    }

    gathered_ = true;
}

const MatrixXd& OLS::X() { gather(); return X_; }
const VectorXd& OLS::y() { gather(); return y_; }

void OLS::solve() {
    if (solved_) return;

    gather();
    const unsigned int n = X_.rows(), k = X_.cols();

    if (n < k) throw RankError("Cannot compute OLS estimates with n < k");

    JacobiSVD<MatrixXd> svd(X_, ComputeThinU | ComputeThinV);
    if (svd.rank() < X_.cols()) throw RankError("Cannot compute OLS estimates: independent data is not full column rank (" + std::to_string(svd.rank()) + " < " + std::to_string(X_.cols()) + ")");

    beta_ = svd.solve(y_);

    residuals_ = y_ - X_ * beta_;

    ssr_ = residuals_.squaredNorm();
    s2_ = X_.cols() == X_.rows() ? std::numeric_limits<double>::quiet_NaN() : ssr_ / (n - k);
    if (model_.hasConstant()) {
        // Centered R^2:
        double sum_y = y_.sum();
        R2_ = 1 - ssr_ / (y_.squaredNorm() - sum_y*sum_y / n);
    }
    else {
        // Uncentered R^2
        R2_ = 1 - ssr_ / y_.squaredNorm();
    }

    DiagonalMatrix<double, Dynamic> S2i(k);
    // We're after: V S^{-2} V', where S is the (diagonal) matrix of SVD singular values and V is
    // the SVD "V" matrix.  Computing the ^{-2} on S is trivial, since it's diagonal: just compute
    // ^{-2} of each element.
    S2i.diagonal() = svd.singularValues().array().square().inverse();
    var_beta_ = s2_ * svd.matrixV() * S2i * svd.matrixV().transpose();

    solved_ = true;
}

const VectorXd& OLS::beta() { solve(); return beta_; }
const MatrixXd& OLS::covariance() { solve(); return var_beta_; }
const VectorXd& OLS::residuals() { solve(); return residuals_; }
const double& OLS::s2() { solve(); return s2_; }
const double& OLS::ssr() { solve(); return ssr_; }

std::ostream& operator<<(std::ostream &os, const OLS &ols) {
    os << "OLS for model: " << ols.model();
    if (ols.solved_) {
        os << "\nResults: beta = " << ols.beta_.transpose() << ", R^2 = " << ols.R2_ << "\n";
        os <<   "   s.e.(beta) = " << ols.var_beta_.diagonal().cwiseSqrt().transpose();
    }
    else {
        os << " (not solved)";
    }
    return os;
}

}}
