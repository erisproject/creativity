#include "creativity/belief/Linear.hpp"
#include <eris/debug.hpp>
#include <eris/Random.hpp>
#include <Eigen/QR>

namespace creativity { namespace belief {

using namespace Eigen;

constexpr double Linear::NONINFORMATIVE_N, Linear::NONINFORMATIVE_S2;

Linear::Linear(
        const Ref<const VectorXd> &beta,
        double s2,
        const Ref<const MatrixXd> &V,
        double n,
        std::shared_ptr<MatrixXd> V_inv,
        std::shared_ptr<MatrixXd> V_chol_L
        )
    : beta_{beta}, s2_{s2}, V_{V}, V_inv_{std::move(V_inv)}, V_chol_L_{std::move(V_chol_L)}, n_{n}, K_(beta_.rows())
{
    // Check that the given matrices conform
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
const MatrixXd& Linear::VcholL() const {
    NO_EMPTY_MODEL;
    if (not V_chol_L_)
        V_chol_L_ = std::make_shared<MatrixXd>(V_.llt().matrixL());
    return *V_chol_L_;
}

const bool& Linear::noninformative() const { NO_EMPTY_MODEL; return noninformative_; }

double Linear::predict(const Ref<const RowVectorXd> &Xi) {
    NO_EMPTY_MODEL;
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

Linear Linear::update(const Ref<const VectorXd> &y, const Ref<const MatrixXd> &X) const {
    NO_EMPTY_MODEL;
    if (y.rows() != X.rows())
        throw std::runtime_error("update(y, X) failed: y and X are non-conformable");
    if (X.cols() != K())
        throw std::runtime_error("update(y, X) failed: X has wrong number of columns");

    MatrixXd Xt = X.transpose();
    MatrixXd XtX = Xt * X;

    // We need the inverse of (V^{-1} + X^\top X) for V_post, but store the value in a shared
    // pointer before inverting because there's a good chance that the next object will need the
    // inverse of the inverse, i.e. by calling its own Vinv()--so, by storing the intermediate value
    // before the inverse, we may be able to avoid inverting an inverse later: so in such a case we
    // save time *and* avoid the numerical precision loss from taking an extra inverse.
    auto V_post_inv = std::make_shared<MatrixXd>(Vinv() + XtX);

    MatrixXd V_post = V_post_inv->colPivHouseholderQr().inverse();

    VectorXd beta_post = V_post * (Vinv() * beta_ + Xt * y);
    double n_post = X.rows();
    if (not noninformative()) n_post += n_;

    VectorXd residualspost = y - X * beta_post;
    VectorXd beta_diff = beta_post - beta_;
    double s2_post = (n_ * s2_ + residualspost.transpose() * residualspost + beta_diff.transpose() * Vinv() * beta_diff) / n_post;

    return Linear(beta_post, s2_post, V_post, n_post, V_post_inv);
}

}}
