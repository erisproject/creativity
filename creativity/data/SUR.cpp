#include "creativity/data/SUR.hpp"
#include "creativity/data/Regression.hpp"
#include "creativity/data/Variable.hpp"
#include <Eigen/SVD>
#include <Eigen/QR>
#include <stdexcept>
#include <string>
#include <ostream>

using namespace Eigen;
namespace creativity { namespace data {

const std::vector<Equation>& SUR::equations() const { return eqs_; }

void SUR::clear() {
    gathered_ = false;
    X_.resize(0, 0);
    y_.resize(0);
    solved_ = false;
    beta_.clear();
    beta_full_.resize(0);
    var_beta_.resize(0, 0);
    residuals_.resize(0);
    ssr_ = s2_ = R2_ = std::numeric_limits<double>::quiet_NaN();
    // FIXME: other vars?
}

void SUR::gather() {
    if (gathered_) return;

    // Since this is an SUR model, each model should have the same number of observations, so check
    // that.  The number of regressors, on the other hand, is the sum of the number of regressors in
    // each equation (which can differ between equations).
    unsigned int per_eq_rows = eqs_.front().depVar().size();
    unsigned int total_cols = 0;
    for (const auto& eq : eqs_) {
        if (eq.depVar().size() != per_eq_rows) throw Variable::SizeError("Invalid SUR: found equations with differing numbers of rows");
        total_cols += eq.numVars();
    }
    unsigned int total_rows = per_eq_rows * eqs_.size();

    y_.resize(total_rows);
    X_.resize(total_rows, total_cols);
    // Second pass: actually populate
    unsigned int r = 0, c = 0;
    for (const auto& eq : eqs_) {
        eq.depVar().populate(y_.segment(r, per_eq_rows));
        for (const Variable &var : eq) {
            var.populate(X_.col(c++).segment(r, per_eq_rows));
        }
        r += per_eq_rows;
    }
    if (r != total_rows or c != total_cols) {
        throw std::runtime_error("Internal error: unexpected number of rows/columns");
    }

    // FIXME: use a sparse matrix?
    
    gathered_ = true;
}

const MatrixXd& SUR::X() { gather(); return X_; }
const VectorXd& SUR::y() { gather(); return y_; }

void SUR::solve() {
    if (solved_) return;

    gather();

    // Step 1: run OLS
    const unsigned int per_eq_rows = eqs_.front().depVar().size(),
          num_eqs = eqs_.size(),
          n = X_.rows(),
          k = X_.cols();

    if (n < k) throw RankError("Cannot compute SUR estimates with n < k");

    JacobiSVD<MatrixXd> svd(X_, ComputeThinU | ComputeThinV);
    if (svd.rank() < X_.cols()) throw RankError("Cannot compute SUR estimates: independent data is not full column rank (" + std::to_string(svd.rank()) + " < " + std::to_string(X_.cols()) + ")");

    // This is the full set of residuals from all the stacked equations:
    MatrixXd res = y_ - X_ * svd.solve(y_);
    // Reshape the matrix so that each column "i" is the residuals from equation i:
    if (res.size() != per_eq_rows * num_eqs) throw std::logic_error("Internal error: residuals has unexpected size");
    res.resize(per_eq_rows, num_eqs);

    // Now get the Sigma estimate from the reshaped residuals, and its inverse
    MatrixXd smallsigmahat = (1.0 / (double) per_eq_rows * res.transpose()) * res;
    MatrixXd smallsigmahatinv = smallsigmahat.fullPivHouseholderQr().inverse();

    // We're aiming at getting:
    // beta = (X^T Z X)^{-1} X^T Z y
    // where Z is sigmahatinv kronecker-product identity; calculate it first:
    MatrixXd bigsigmahat = MatrixXd::Zero(per_eq_rows * num_eqs, per_eq_rows * num_eqs);
    for (unsigned int i = 0; i < num_eqs; i++) {
        for (unsigned int j = 0; j < num_eqs; j++) {
            bigsigmahat.block(i*per_eq_rows, j*per_eq_rows, per_eq_rows, per_eq_rows).diagonal().setConstant(smallsigmahatinv(i,j));
        }
    }

    // Precalculate the X^T Z part, since we need it twice:
    MatrixXd xtsig = X_.transpose() * bigsigmahat;

    beta_full_ = (xtsig * X_).fullPivHouseholderQr().solve(xtsig * y_);

    for (unsigned int i = 0, k = 0; i < num_eqs; k += eqs_[i++].numVars()) {
        beta_.emplace_back(beta_full_, k, eqs_[i].numVars());
    }

/*   FIXME:
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
*/

    solved_ = true;
}

const std::vector<VectorBlock<VectorXd>>& SUR::beta() { solve(); return beta_; }
const MatrixXd& SUR::covariance() { solve(); return var_beta_; }
const VectorXd& SUR::residuals() { solve(); return residuals_; }
const double& SUR::s2() { solve(); return s2_; }
const double& SUR::ssr() { solve(); return ssr_; }

std::ostream& operator<<(std::ostream &os, const SUR &sur) {
    os << "SUR system:\n";
    for (auto &eq : sur.equations()) {
        os << "\t" << eq << "\n";
    }
    if (sur.solved_) {
        os << "\nBeta:\n";
        for (auto &beta : sur.beta_) {
            os << "\t" << beta.transpose() << "\n";
        }
    }
    else {
        os << " (not solved)";
    }
    return os;
}

}}
