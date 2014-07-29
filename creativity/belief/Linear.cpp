#include "creativity/belief/Linear.hpp"
#include <eris/debug.hpp>
#include <Eigen/QR>

namespace creativity { namespace belief {

using namespace Eigen;

constexpr double Linear::NONINFORMATIVE_N, Linear::NONINFORMATIVE_S2;

Linear::Linear(
        const Ref<const VectorXd> &beta,
        double s2,
        const Ref<const MatrixXd> &V,
        double n,
        std::shared_ptr<MatrixXd> V_inv
        )
    : beta_{beta}, s2_{s2}, V_{V}, V_inv_{std::move(V_inv)}, n_{n}, K_(beta_.rows())
{
    // Check that the given matrices conform
    auto k = K();
    if (k < 1) throw std::logic_error("Linear model requires at least one parameter");
    if (V_.rows() != V_.cols()) throw std::logic_error("Linear requires square V matrix");
    if (k != V_.rows()) throw std::logic_error("Linear requires beta and V of same number of rows");
    if (V_inv_ and (V_inv_->rows() != V_inv_->cols() or V_inv_->rows() != k))
        throw std::logic_error("Linear constructed with invalid V_inv");
    auto fixed = fixedModelSize();
    if (fixed and k != fixed) throw std::logic_error("Linear model constructed with incorrect number of model parameters");
}

Linear::Linear(unsigned int K) :
    beta_{VectorXd::Zero(K)}, s2_{NONINFORMATIVE_S2}, V_{MatrixXd::Identity(K, K)},
    n_{NONINFORMATIVE_N}, noninformative_{true}, K_{K}
{
    if (K < 1) throw std::logic_error("Linear model requires at least one parameter");
    auto fixed = fixedModelSize();
    if (fixed and K != fixed) throw std::logic_error("Linear model constructed with incorrect number of model parameters");
}

unsigned int Linear::fixedModelSize() const { return 0; }

#define NO_EMPTY_MODEL if (K_ == 0) { throw std::logic_error("Cannot use default constructed model object as a model"); }

const VectorXd& Linear::beta() const { NO_EMPTY_MODEL; return beta_; }
const double& Linear::s2() const { NO_EMPTY_MODEL; return s2_; }
const double& Linear::n() const { NO_EMPTY_MODEL; return n_; }
const MatrixXd& Linear::V() const { NO_EMPTY_MODEL; return V_; }
const MatrixXd& Linear::Vinv() const {
    NO_EMPTY_MODEL;
    if (not V_inv_)
        V_inv_ = std::make_shared<MatrixXd>(V_.colPivHouseholderQr().inverse());
    return *V_inv_;
}

const bool& Linear::noninformative() const { NO_EMPTY_MODEL; return noninformative_; }

double Linear::predict(const Ref<const RowVectorXd> &Xi) const {
    NO_EMPTY_MODEL;
    // FIXME: this is wrong; prediction should be using some sort of sampling
    return Xi * beta_;
}

std::ostream& operator<<(std::ostream &os, const Linear &b) {
    if (b.K() == 0)
        os << "Linear model with no parameters (default constructed)";
    else
        os << "Linear model with " << b.K() << " parameters, beta_ =\n" << b.beta_;
    return os;
}

void Linear::verifyParameters() const { NO_EMPTY_MODEL; }

Linear Linear::update(const Ref<const VectorXd> &y, const Ref<const MatrixXd> &X) const {
    NO_EMPTY_MODEL;
    if (y.rows() != X.rows())
        throw std::runtime_error("update(y, X) failed: y and X are non-conformable");
    if (X.cols() != K())
        throw std::runtime_error("update(y, X) failed: X has wrong number of columns");

    MatrixXd Xt = X.transpose();
    MatrixXd XtX = Xt * X;

    // Store (V^{-1} + X^\top X) in a shared pointer so that we can pass it along as V^{-1}
    // to the new Linear object so that if *it* needs an inverse, it can just use this
    // instead of having to reinvert to get back to this value.  This saves some
    // computational time and, more importantly, gives more accurate numerical results.
    auto V_post_inv = std::make_shared<MatrixXd>(Vinv() + XtX);

    MatrixXd V_post = V_post_inv->colPivHouseholderQr().inverse();

    VectorXd beta_post = V_post * (Vinv() * beta_ + Xt * y);
    double n_post = X.rows();
    if (not noninformative()) n_post += n_;

    VectorXd residuals = y - X * beta_;
    double s2_post = (n_ * s2_ + residuals.transpose() * (X * V_ * Xt + MatrixXd::Identity(X.rows(), X.rows())) * residuals) / n_post;

    return Linear(beta_post, s2_post, V_post, n_post, V_post_inv);
}

}}
