#include "creativity/belief/Linear.hpp"
#include <eris/debug.hpp>
#include <eris/Random.hpp>
#include <Eigen/QR>

namespace creativity { namespace belief {

using namespace Eigen;

constexpr double Linear::NONINFORMATIVE_N, Linear::NONINFORMATIVE_S2;

Linear::Linear(
        const Ref<const VectorXd> beta,
        double s2,
        const Ref<const MatrixXd> V,
        double n,
        std::shared_ptr<MatrixXd> V_inv,
        std::shared_ptr<MatrixXd> V_chol_L
        )
    : beta_(beta), s2_{s2}, V_(V), V_inv_{std::move(V_inv)}, V_chol_L_{std::move(V_chol_L)}, n_{n}, K_(beta_.rows())
{
    // Check that the given matrices conform
    checkLogic();
}

Linear::Linear(
        const std::vector<double> &beta,
        double s2,
        const std::vector<double> &V,
        double n
        )
    : beta_(beta.size()), s2_(s2), V_(beta.size(), beta.size()), n_(n), K_(beta.size())
{
    checkLogic();
    if (V.size() != K_*(K_+1)/2)
        throw std::logic_error("Linear vector constructor called with invalid V vector for model with K=" + std::to_string(K_) +
                " (expected " + std::to_string(K_*(K_+1)/2) + " lower triangle elements, received " + std::to_string(V.size()) + ")");

    for (unsigned int k = 0; k < K_; k++) beta_[k] = beta[k];
    for (unsigned int r = 0, i = 0; r < K_; r++) for (unsigned int c = 0; c <= r; c++, i++) {
        const double &vij = V[i];
        V_(r,c) = vij;
        if (r != c) V_(c,r) = vij;
    }
}

Linear::Linear(unsigned int K) :
    beta_{VectorXd::Zero(K)}, s2_{NONINFORMATIVE_S2}, V_{MatrixXd::Identity(K, K)},
    n_{NONINFORMATIVE_N}, noninformative_{true}, K_{K}
{
    if (K < 1) throw std::logic_error("Linear model requires at least one parameter");
    auto fixed = fixedModelSize();
    if (fixed and K != fixed) throw std::logic_error("Linear model constructed with incorrect number of model parameters");
}


void Linear::checkLogic() {
    auto k = K();
    if (k < 1) throw std::logic_error("Linear model requires at least one parameter");
    if (V_.rows() != V_.cols()) throw std::logic_error("Linear requires square V matrix");
    if (k != V_.rows()) throw std::logic_error("Linear requires beta and V of same number of rows");
    if (V_inv_ and (V_inv_->rows() != V_inv_->cols() or V_inv_->rows() != k))
        throw std::logic_error("Linear constructed with invalid V_inv");
    if (V_chol_L_ and (V_chol_L_->rows() != V_chol_L_->cols() or V_chol_L_->rows() != k))
        throw std::logic_error("Linear constructed with invalid V_chol_L");
    auto fixed = fixedModelSize();
    if (fixed and k != fixed) throw std::logic_error("Linear model constructed with incorrect number of model parameters");
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
const MatrixXd& Linear::VcholL() const {
    NO_EMPTY_MODEL;
    if (not V_chol_L_)
        V_chol_L_ = std::make_shared<MatrixXd>(V_.llt().matrixL());
    return *V_chol_L_;
}

const bool& Linear::noninformative() const { NO_EMPTY_MODEL; return noninformative_; }

double Linear::predict(const Ref<const RowVectorXd> &Xi) {
    NO_EMPTY_MODEL;
    if (noninformative_)
        throw std::logic_error("Cannot call predict() on noninformative model");
    return Xi * beta_;
}

const VectorXd& Linear::draw() {
    NO_EMPTY_MODEL;

    if (last_draw_.size() != K_ + 1) last_draw_.resize(K_ + 1);

    auto &rng = eris::Random::rng();

    // beta is distributed as t(beta, s^2*V, n)
    
    // That can be generated as beta + y*sqrt(n/q) where y ~ N(0, s^2*V), and q ~ chisq(n)

    // To generate y ~ N(0, SIGMA), generate N(0,1) and multiple by L from the cholesky
    // decomposition of SIGMA, or in other words, s*L where LL'=V (so SIGMA=s^2*V)
    VectorXd y(K_);
    std::normal_distribution<double> stdnorm(0, 1);
    for (unsigned int i = 0; i < K_; i++) y[i] = stdnorm(rng);

    std::chi_squared_distribution<double> rchisqn(n_);

    last_draw_.head(K_) = beta_ + sqrt(s2_ * n_ / rchisqn(rng)) * VcholL() * y;

    // h has distribution Gamma(n/2, 2/(n*s^2)), and s^2 is the inverse of h:
    last_draw_[K_] = 1.0 / (std::gamma_distribution<double>(n_/2, 2/(s2_*n_))(rng));

    return last_draw_;
}

const VectorXd& Linear::lastDraw() const {
    return last_draw_;
}

void Linear::discard(unsigned int) {
    NO_EMPTY_MODEL;
    // Default implementation is a no-op since the default draw() distributions are already
    // independent.
}

void Linear::discardForce(unsigned int burn) {
    NO_EMPTY_MODEL;
    for (unsigned int i = 0; i < burn; i++)
        draw();
}

std::ostream& operator<<(std::ostream &os, const Linear &b) {
    if (b.K() == 0)
        os << "Linear model with no parameters (default constructed)";
    else
        os << "Linear model with " << b.K() << " parameters, beta_ =\n" << b.beta_;
    return os;
}

void Linear::verifyParameters() const { NO_EMPTY_MODEL; }

// Called on an lvalue object, creates a new object with *this as prior
Linear Linear::update(const Ref<const VectorXd> &y, const Ref<const MatrixXd> &X) const & {
    Linear updated(*this);
    updated.updateInPlace(y, X);
    return updated;
}

// Called on rvalue, so just update *this as needed, using itself as the prior, then return std::move(*this)
Linear Linear::update(const Ref<const VectorXd> &y, const Ref<const MatrixXd> &X) && {
    updateInPlace(y, X);
    return std::move(*this);
}

void Linear::updateInPlace(const Ref<const VectorXd> &y, const Ref<const MatrixXd> &X) {
    NO_EMPTY_MODEL;
    if (X.cols() != K())
        throw std::logic_error("update(y, X) failed: X has wrong number of columns");
    if (y.rows() != X.rows())
        throw std::logic_error("update(y, X) failed: y and X are non-conformable");

    reset();

    if (y.rows() == 0) // Nothing to update!
        return;

    MatrixXd Xt = X.transpose();
    MatrixXd XtX = Xt * X;

    // We need the inverse of (V^{-1} + X^\top X) for V_post, but store the value in a shared
    // pointer before inverting because there's a good chance that the next object will need the
    // inverse of the inverse, i.e. by calling its own Vinv()--so, by storing the intermediate value
    // before the inverse, we may be able to avoid inverting an inverse later: so in such a case we
    // save time *and* avoid the numerical precision loss from taking an extra inverse.
    auto V_post_inv = std::make_shared<MatrixXd>(Vinv() + XtX);
    V_ = V_post_inv->colPivHouseholderQr().inverse();
    VectorXd beta_post = V_ * (Vinv() * beta_ + Xt * y);

    double n_prior = noninformative() ? 0 : n_;
    n_ = n_prior + X.rows();

    VectorXd residualspost = y - X * beta_post;
    VectorXd beta_diff = beta_post - beta_;
    beta_ = std::move(beta_post);
    s2_ = (n_prior * s2_ + residualspost.squaredNorm() + beta_diff.transpose() * Vinv() * beta_diff) / n_;

    V_inv_ = std::move(V_post_inv);
    beta_ = std::move(beta_post);
    if (V_chol_L_) V_chol_L_.reset(); // This will have to be recalculated
    if (noninformative_) noninformative_ = false; // If we just updated a noninformative model, we aren't noninformative anymore
}

Linear Linear::weaken(const double precision_scale) const & {
    Linear weakened(*this);
    weakened.weakenInPlace(precision_scale);
    return weakened;
}

Linear Linear::weaken(const double precision_scale) && {
    weakenInPlace(precision_scale);
    return std::move(*this);
}

void Linear::weakenInPlace(const double precision_scale) {
    if (precision_scale <= 0 or precision_scale > 1)
        throw std::logic_error("weaken() called with invalid precision multiplier (not in (0,1])");

    reset();

    if (noninformative() or precision_scale == 1.0) // Nothing to do here
        return;

    // If V_inv_ is set, scale it (because this is much cheaper than inverting the scaled V_ later)
    if (V_inv_) {
        // If we're the unique owner of the inverse matrix, scale it directly
        if (V_inv_.unique())
            *V_inv_ *= precision_scale;
        // Otherwise we have to make a copy (because we don't want to change the original!)
        else
            V_inv_ = std::make_shared<Eigen::MatrixXd>(*V_inv_ * precision_scale);
    }

    // Likewise for the Cholesky decomposition
    if (V_chol_L_) {
        if (V_chol_L_.unique())
            *V_chol_L_ /= std::sqrt(precision_scale);
        else
            V_chol_L_ = std::make_shared<Eigen::MatrixXd>(*V_chol_L_ / std::sqrt(precision_scale));
    }

    V_ /= precision_scale;

    return;
}

void Linear::reset() {
    if (last_draw_.size() > 0) last_draw_.resize(0);
}

}}
