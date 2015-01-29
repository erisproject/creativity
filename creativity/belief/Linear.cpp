#include "creativity/belief/Linear.hpp"
#include <eris/Random.hpp>
#include <Eigen/QR>

namespace creativity { namespace belief {

using namespace Eigen;
using eris::Random;

constexpr double Linear::NONINFORMATIVE_N, Linear::NONINFORMATIVE_S2, Linear::NONINFORMATIVE_Vc;

Linear::Linear(
        const Ref<const VectorXd> beta,
        double s2,
        const Ref<const MatrixXd> V,
        double n)
    : beta_(beta), s2_{s2}, V_(V.selfadjointView<Lower>()), n_{n}, K_(beta_.rows())
{
    // Check that the given matrices conform
    checkLogic();
}

Linear::Linear(unsigned int K, const Ref<const MatrixXdR> noninf_X, const Ref<const VectorXd> noninf_y) :
    beta_{VectorXd::Zero(K)}, s2_{NONINFORMATIVE_S2}, V_{NONINFORMATIVE_Vc * MatrixXd::Identity(K, K)},
    n_{NONINFORMATIVE_N}, noninformative_{true}, K_{K}
{
    if (K < 1) throw std::logic_error("Linear model requires at least one parameter");
    auto fixed = fixedModelSize();
    if (fixed and K != fixed) throw std::logic_error("Linear model constructed with incorrect number of model parameters");

    if (noninf_X.rows() != noninf_y.rows()) throw std::logic_error("Partially informed model construction error: X.rows() != y.rows()");
    if (noninf_X.rows() > 0) {
        if (noninf_X.cols() != K_) throw std::logic_error("Partially informed model construction error: X.cols() != K");
        noninf_X_.reset(new MatrixXdR(noninf_X));
        noninf_y_.reset(new VectorXd(noninf_y));
    }
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
        V_inv_ = std::make_shared<MatrixXd>(V_.fullPivHouseholderQr().inverse().selfadjointView<Lower>());
    return *V_inv_;
}
const MatrixXd& Linear::VcholL() const {
    NO_EMPTY_MODEL;
    if (not V_chol_L_)
        V_chol_L_ = std::make_shared<MatrixXd>(V_.selfadjointView<Lower>().llt().matrixL());
    return *V_chol_L_;
}
const MatrixXd& Linear::VcholLinv() const {
    NO_EMPTY_MODEL;
    if (not V_chol_L_inv_)
        V_chol_L_inv_ = std::make_shared<MatrixXd>(VcholL().fullPivHouseholderQr().inverse());
    return *V_chol_L_inv_;
}

const bool& Linear::noninformative() const { NO_EMPTY_MODEL; return noninformative_; }

const std::vector<std::string>& Linear::names() const {
    if (not beta_names_ or beta_names_->size() != K()) {
        beta_names_.reset(new std::vector<std::string>);
        for (unsigned int i = 0; i < K_; i++) {
            beta_names_->push_back(std::to_string(i));
        }
    }
    return *beta_names_;
}

void Linear::names(const std::vector<std::string> &names) {
    if (names.empty()) {
        beta_names_.reset();
        beta_names_default_ = true;
        return;
    }
    else if (names.size() != K_) throw std::domain_error("Linear::names(): given names vector is not of size K");

    beta_names_.reset(new std::vector<std::string>());
    beta_names_->reserve(K_);
    for (const auto &n : names) beta_names_->push_back(n);
    beta_names_default_ = false;
}

double Linear::predict(const Ref<const RowVectorXd> &Xi) {
    NO_EMPTY_MODEL;
    if (noninformative_)
        throw std::logic_error("Cannot call predict() on noninformative model");
    return Xi * beta_;
}

const VectorXd& Linear::draw() {
    NO_EMPTY_MODEL;

    if (last_draw_.size() != K_ + 1) last_draw_.resize(K_ + 1);

    auto &rng = Random::rng();

    // (beta,h) is distributed as a normal-gamma(beta, V, s2^{-1}, n), in Koop's Gamma distribution
    // notation, or NG(beta, V, n/2, 2*s2^{-1}/n) in the more common G(shape,scale) notation
    // (which std::gamma_distribution wants).
    //
    // Proof:
    // Let $G_{k\theta}(k,\theta)$ be the shape ($k$), scale ($\theta$) notation.  This has mean $k\theta$ and
    // variance $k\theta^2$.
    //
    // Let $G_{Koop}(\mu,\nu)$ be Koop's notation, where $\mu$ is the mean and $\nu$ is the degrees of
    // freedom, which has variance $\frac{2\mu^2}{\nu}$.  Equating means and variances:
    //
    // \[
    //     k\theta = \mu
    //     k\theta^2 = \frac{2\mu^2}{\nu}
    //     \theta = \frac{2\mu}{\nu}
    //     k = \frac{2}{\nu}
    // \]
    // where the third equation follows from the first divided by the second, and fourth follows
    // from the first divided by the third.  Thus
    // \[
    //     G_{Koop}(\mu,\nu) = G_{k\theta}(\frac{2}{\nu},\frac{2\mu}{\nu})
    // \]

    // To draw this, first draw a gamma-distributed "h" value (store its inverse)
    last_draw_[K_] = 1.0 / std::gamma_distribution<double>(n_/2, 2/(s2_*n_))(rng);

    // Now use that to draw a multivariate normal conditional on h, with mean beta and variance
    // h^{-1} V; this is the beta portion of the draw:
    last_draw_.head(K_) = multivariateNormal(beta_, VcholL(), std::sqrt(last_draw_[K_]));

    return last_draw_;
}

VectorXd Linear::multivariateNormal(const Ref<const VectorXd> &mu, const Ref<const MatrixXd> &L, double s) {
    if (mu.rows() != L.rows() or L.rows() != L.cols())
        throw std::logic_error("multivariateNormal() called with non-conforming mu and L");

    // To draw such a normal, we need the lower-triangle Cholesky decomposition L of V, and a vector
    // of K random \f$N(\mu=0, \sigma^2=h^{-1})\f$ values.  Then \f$beta + Lz\f$ yields a \f$beta\f$
    // draw of the desired distribution.
    VectorXd z(mu.size());
    for (unsigned int i = 0; i < z.size(); i++) z[i] = Random::rstdnorm();

    return mu + L * (s * z);
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
    return os << (std::string) b;
}

Linear::operator std::string() const {
    std::ostringstream summary;
    summary << display_name();
    if (K_ == 0) summary << " model with no parameters (default constructed)";
    else {
        if (noninformative()) summary << " (noninformative)";
        summary << " model: K=" << K_ << ", n=" << n_ << ", s2=" << s2_;
        if (not beta_names_default_) {
            summary << "\n  X cols:";
            for (auto &n : names()) {
                if (n.find_first_of(" \t\n") != n.npos) summary << " {" << n << "}";
                else summary << " " << n;
            }
        }
        summary <<
            "\n  beta = " << beta_.transpose().format(IOFormat(StreamPrecision, 0, ", ")) <<
            "\n  V = " << V_.format(IOFormat(6, 0, " ", "\n      ")) << "\n";
    }
    return summary.str();
}

std::string Linear::display_name() const { return "Linear"; }

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
        throw std::logic_error("update(y, X) failed: X has wrong number of columns (expected " + std::to_string(K()) + ", got " + std::to_string(X.cols()) + ")");
    if (y.rows() != X.rows())
        throw std::logic_error("update(y, X) failed: y and X are non-conformable");

    reset();

    if (y.rows() == 0) // Nothing to update!
        return;

    if (noninformative_) {
        if (not noninf_X_ or noninf_X_->rows() == 0) noninf_X_.reset(new MatrixXdR(X.rows(), K()));
        else if (not noninf_X_.unique()) {
            auto old_x = noninf_X_;
            noninf_X_.reset(new MatrixXdR(old_x->rows() + X.rows(), K()));
            noninf_X_->topRows(old_x->rows()) = *old_x;
        }
        else noninf_X_->conservativeResize(noninf_X_->rows() + X.rows(), K());

        if (not noninf_y_ or noninf_y_->rows() == 0) noninf_y_.reset(new VectorXd(y.rows()));
        else if (not noninf_y_.unique()) {
            auto old_y = noninf_y_;
            noninf_y_.reset(new VectorXd(old_y->rows() + y.rows()));
            noninf_y_->head(old_y->rows()) = *old_y;
        }
        else noninf_y_->conservativeResize(noninf_y_->rows() + y.rows());

        noninf_X_->bottomRows(X.rows()) = X;
        noninf_y_->tail(y.rows()) = y;

        if (noninf_X_->rows() > K()) {
            MatrixXd XtX = noninf_X_->transpose() * *noninf_X_;
            auto qr = XtX.fullPivHouseholderQr();
            if (qr.rank() >= K()) {
                V_ = qr.inverse().selfadjointView<Lower>();
                V_inv_.reset(new MatrixXd(
#ifdef EIGEN_HAVE_RVALUE_REFERENCES
                        std::move
#endif
                        (XtX)));
                beta_ = V_ * (noninf_X_->transpose() * *noninf_y_);
                n_ = X.rows();
                s2_ = (*noninf_y_ - *noninf_X_ * beta_).squaredNorm() / (n_ - K());
    
                // These are unlikely to be set, but just in case:
                if (V_chol_L_) V_chol_L_.reset();
                if (V_chol_L_inv_) V_chol_L_inv_.reset(); 

                noninf_X_.reset();
                noninf_y_.reset();
                noninformative_ = false; // We aren't noninformative anymore!
            }
        }
    }
    else {
        // Otherwise we were already informative, so just pass the data along.
        updateInPlaceInformative(y, X);
    }
}

void Linear::updateInPlaceInformative(const Ref<const VectorXd> &y, const Ref<const MatrixXd> &X) {
    MatrixXd Xt = X.transpose();
    MatrixXd XtX = Xt * X;

    // We need the inverse of (V^{-1} + X^\top X) for V_post, but store the value in a shared
    // pointer before inverting because there's a good chance that the next object will need the
    // inverse of the inverse, i.e. by calling its own Vinv()--so, by storing the intermediate value
    // before the inverse, we may be able to avoid inverting an inverse later: so in such a case we
    // save time *and* avoid the numerical precision loss from taking an extra inverse.
    auto V_post_inv = std::make_shared<MatrixXd>(Vinv() + XtX);
    V_ = V_post_inv->fullPivHouseholderQr().inverse();
    VectorXd beta_post(V_ * (Vinv() * beta_ + Xt * y));

    double n_prior = n_;
    n_ += X.rows();

    VectorXd residualspost = y - X * beta_post;
    VectorXd beta_diff = beta_post - beta_;
    s2_ = (n_prior * s2_ + residualspost.squaredNorm() + beta_diff.transpose() * Vinv() * beta_diff) / n_;
    beta_ =
#ifdef EIGEN_HAVE_RVALUE_REFERENCES
        std::move(beta_post);
#else
        beta_post;
#endif

    V_inv_ = std::move(V_post_inv);
    // The decompositions will have to be recalculated:
    if (V_chol_L_) V_chol_L_.reset();
    if (V_chol_L_inv_) V_chol_L_inv_.reset(); 
}

Linear Linear::weaken(const double stdev_scale) const & {
    Linear weakened(*this);
    weakened.weakenInPlace(stdev_scale);
    return weakened;
}

Linear Linear::weaken(const double stdev_scale) && {
    weakenInPlace(stdev_scale);
    return std::move(*this);
}

void Linear::weakenInPlace(const double stdev_scale) {
    if (stdev_scale < 1)
        throw std::logic_error("weaken() called with invalid stdev multiplier " + std::to_string(stdev_scale) + " < 1");

    reset();

    if (noninformative() or stdev_scale == 1.0) // Nothing to do here
        return;

    if (noninf_X_) {
        // Partially informed model
        if (not noninf_X_.unique()) noninf_X_.reset(new MatrixXdR(*noninf_X_ / stdev_scale));
        else *noninf_X_ /= stdev_scale;
        if (not noninf_y_.unique()) noninf_y_.reset(new VectorXd(*noninf_y_ / stdev_scale));
        else *noninf_y_ /= stdev_scale;
        return;
    }

    const double var_scale = stdev_scale*stdev_scale;

    // If V_inv_ is set, scale it (because this is much cheaper than inverting the scaled V_ later)
    if (V_inv_) {
        // If we're the unique owner of the inverse matrix, scale it directly
        if (V_inv_.unique())
            *V_inv_ /= var_scale;
        // Otherwise we have to make a copy (because we don't want to change the original!)
        else
            V_inv_ = std::make_shared<Eigen::MatrixXd>(*V_inv_ / var_scale);
    }

    // Likewise for the Cholesky decomposition (and its inverse)
    if (V_chol_L_) {
        if (V_chol_L_.unique())
            *V_chol_L_ *= stdev_scale;
        else
            V_chol_L_.reset(new Eigen::MatrixXd(*V_chol_L_ * stdev_scale));
    }
    if (V_chol_L_inv_) {
        if (V_chol_L_inv_.unique())
            *V_chol_L_inv_ /= stdev_scale;
        else
            V_chol_L_inv_.reset(new Eigen::MatrixXd(*V_chol_L_inv_ / stdev_scale));
    }

    // And of course V gets scaled
    V_ *= var_scale;

    return;
}

const MatrixXdR& Linear::noninfXData() const {
    if (not noninformative_) throw std::logic_error("noninfXData() cannot be called on a fully-informed model");

    if (not noninf_X_) const_cast<std::shared_ptr<MatrixXdR>&>(noninf_X_).reset(new MatrixXdR(0, K()));

    return *noninf_X_;
}

const VectorXd& Linear::noninfYData() const {
    if (not noninformative_) throw std::logic_error("noninfYData() cannot be called on a fully-informed model");

    if (not noninf_y_) const_cast<std::shared_ptr<VectorXd>&>(noninf_y_).reset(new VectorXd(0));

    return *noninf_y_;
}

void Linear::reset() {
    if (last_draw_.size() > 0) last_draw_.resize(0);
}

Linear::draw_failure::draw_failure(const std::string &what) : std::runtime_error(what) {}

Linear::draw_failure::draw_failure(const std::string &what, const Linear &model) : std::runtime_error(what + "\n" + (std::string) model) {}

}}
