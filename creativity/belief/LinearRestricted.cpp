// Eigen uses deprecated binder1st and binder2nd; until that is fixed, ignore the generated warning:
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include "creativity/belief/LinearRestricted.hpp"
#include <cmath>
#include <Eigen/QR>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>

#include <eris/debug.hpp>

using namespace Eigen;
using eris::Random;

namespace creativity { namespace belief {

LinearRestricted::RestrictionProxy LinearRestricted::lowerBound(size_t k) {
    return RestrictionProxy(*this, k, false);
}
const LinearRestricted::RestrictionProxy LinearRestricted::lowerBound(size_t k) const {
    return RestrictionProxy(const_cast<LinearRestricted&>(*this), k, false);
}
LinearRestricted::RestrictionProxy LinearRestricted::upperBound(size_t k) {
    return RestrictionProxy(*this, k, true);
}
const LinearRestricted::RestrictionProxy LinearRestricted::upperBound(size_t k) const {
    return RestrictionProxy(const_cast<LinearRestricted&>(*this), k, true);
}
LinearRestricted::RestrictionIneqProxy LinearRestricted::restrict(size_t k) {
    return RestrictionIneqProxy(*this, k);
}
const LinearRestricted::RestrictionIneqProxy LinearRestricted::restrict(size_t k) const {
    return RestrictionIneqProxy(const_cast<LinearRestricted&>(*this), k);
}

void LinearRestricted::allocateRestrictions(size_t more) {
    // Increment by K_ rows at a time.  (This is fairly arbitrary, but at least the number of
    // restrictions will typically be correlated with the number of regressors)
    size_t rows = K_ * (size_t) std::ceil((restrict_size_ + more) / (double) K_);
    if (restrict_select_.cols() != K_) restrict_select_.conservativeResize(rows, K_);
    else if ((size_t) restrict_select_.rows() < restrict_size_ + more) {
        restrict_select_.conservativeResize(rows, NoChange);
    }
    if ((size_t) restrict_values_.size() < restrict_size_ + more) {
        restrict_values_.conservativeResize(rows);
    }
}

void LinearRestricted::addRestriction(const Ref<const RowVectorXd> &R, double r) {
    if (R.size() != K_) throw std::logic_error("Unable to add linear restriction: R does not have size K");
    allocateRestrictions(1);
    restrict_values_[restrict_size_] = r;
    restrict_select_.row(restrict_size_) = R;
    restrict_size_++;
    reset();
}

void LinearRestricted::addRestrictionGE(const Ref<const RowVectorXd> &R, double r) {
    return addRestriction(-R, -r);
}

void LinearRestricted::addRestrictions(const Ref<const MatrixXd> &R, const Ref<const VectorXd> &r) {
    if (R.cols() != K_) throw std::logic_error("Unable to add linear restrictions: R does not have K columns");
    auto num_restr = R.rows();
    if (num_restr != r.size()) throw std::logic_error("Unable to add linear restrictions: different number of rows in R and r");
    allocateRestrictions(num_restr);
    restrict_values_.segment(restrict_size_, num_restr) = r;
    restrict_select_.middleRows(restrict_size_, num_restr) = R;
    restrict_size_ += num_restr;
    reset();
}

void LinearRestricted::addRestrictionsGE(const Ref<const MatrixXd> &R, const Ref<const VectorXd> &r) {
    return addRestrictions(-R, -r);
}

void LinearRestricted::clearRestrictions() {
    restrict_size_ = 0;
    reset();
    // Don't need to resize restrict_*, the values aren't used if restrict_size_ = 0
}

void LinearRestricted::discard(unsigned int burn) {
    Linear::discard(burn); // Pass it up (mainly for the no-empty-model check)
    mean_beta_draws_ = 0; // Reset the number of beta draws so that the next predict() redraws
}

void LinearRestricted::discardForce(unsigned int burn) {
    Linear::discardForce(burn);
    mean_beta_draws_ = 0;
}

void LinearRestricted::reset() {
    Linear::reset();
    mean_beta_draws_ = 0;
    draw_rejection_discards_last = 0;
    draw_rejection_discards = 0;
    draw_rejection_success = 0;
    gibbs_D_.reset();
    gibbs_last_z_.reset();
    gibbs_last_sigma_ = std::numeric_limits<double>::signaling_NaN();
    gibbs_draws_ = 0;
}

Block<const Matrix<double, Dynamic, Dynamic, RowMajor>, Dynamic, Dynamic, true> LinearRestricted::R() const {
    return restrict_select_.topRows(restrict_size_);
}

VectorBlock<const VectorXd> LinearRestricted::r() const {
    return restrict_values_.head(restrict_size_);
}

const VectorXd& LinearRestricted::draw() {
    return draw(draw_mode);
}

const VectorXd& LinearRestricted::draw(DrawMode m) {
    // If they explicitly want rejection draw, do it:
    if (m == DrawMode::Rejection)
        return drawRejection();

    // In auto mode we might try rejection sampling first
    if (m == DrawMode::Auto) { // Need success rate < 0.1 to switch, with at least 20 samples.
        long draw_rej_samples = draw_rejection_success + draw_rejection_discards;
        double success_rate = draw_rejection_success / (double)draw_rej_samples;
        if (success_rate < draw_auto_min_success_rate and draw_rej_samples >= draw_rejection_max_discards) {
            // Too many failures, switch to Gibbs sampling
            m = DrawMode::Gibbs;
        }
        else {
            // Figure out how many sequential failed draws we'd need to hit the failure threshold,
            // then try up to that many times.
            long draw_tries = std::max<long>(
                    std::ceil(draw_rejection_success / draw_auto_min_success_rate),
                    draw_rejection_max_discards)
                - draw_rej_samples;
            try {
                return drawRejection(draw_tries);
            }
            catch (draw_failure&) {
                // Draw failure; switch to Gibbs
                m = DrawMode::Gibbs;
            }
        }
    }

    // Either Gibbs was requested explicitly, or Auto was asked for and rejection sampling failed
    // (either we tried it just now, or previous draws have too low a success rate), so use Gibbs
    return drawGibbs();
}

void LinearRestricted::gibbsInitialize(const Ref<const VectorXd> &initial, unsigned long max_tries) {
    constexpr double overshoot = 1.5;

    if (initial.size() < K_ or initial.size() > K_+1) throw std::logic_error("LinearRestricted::gibbsInitialize() called with invalid initial vector (initial.size() != K())");

    gibbs_draws_ = 0;

    const MatrixXd &A = VcholLinv();
    auto &rng = Random::rng();

    if (restrict_size_ == 0) {
        // No restrictions, nothing special to do!
        if (not gibbs_last_z_ or gibbs_last_z_->size() != K_) gibbs_last_z_.reset(new VectorXd(K_));
        *gibbs_last_z_ = A / std::sqrt(s2_) * (initial.head(K_) - beta_);
    }
    else {
        // Otherwise we'll start at the initial value and update
        VectorXd adjusted(initial.head(K_));
        std::vector<size_t> violations;
        violations.reserve(restrict_size_);

        for (unsigned long trial = 1; ; trial++) {
            violations.clear();
            VectorXd v = restrict_select_.topRows(restrict_size_) * adjusted - restrict_values_.head(restrict_size_);
            for (size_t i = 0; i < restrict_size_; i++) {
                if (v[i] > 0) violations.push_back(i);
            }
            if (violations.size() == 0) break; // Restrictions satisfied!
            else if (trial >= max_tries) { // Too many tries: give up
                // Clear gibbs_last_ on failure (so that the next drawGibbs() won't try to use it)
                gibbs_last_z_.reset();
                throw constraint_failure("gibbsInitialize() couldn't find a way to satisfy the model constraints");
            }

            // Select a random constraint to fix:
            size_t fix = violations[std::uniform_int_distribution<size_t>(0, violations.size()-1)(rng)];

            // Aim at the nearest point on the boundary and overshoot (by 50%):
            adjusted += -overshoot * v[fix] / restrict_select_.row(fix).squaredNorm() * restrict_select_.row(fix).transpose();
        }

        if (not gibbs_last_z_ or gibbs_last_z_->size() != K_) gibbs_last_z_.reset(new VectorXd(K_));
        *gibbs_last_z_ = A / std::sqrt(s2_) * (adjusted - beta_);
    }

    // Don't set the last sigma draw; drawGibbs will do that.
    gibbs_last_sigma_ = std::numeric_limits<double>::signaling_NaN();
}

const VectorXd& LinearRestricted::drawGibbs() {
    last_draw_mode = DrawMode::Gibbs;

    const MatrixXd &Ainv = VcholL();

    if (not gibbs_D_) gibbs_D_.reset(
            restrict_size_ == 0
            ? new Matrix<double, Dynamic, Dynamic, RowMajor>(0, K_) // It's rather pointless to use drawGibbs() with no restrictions, but allow it for debugging purposes
            : new Matrix<double, Dynamic, Dynamic, RowMajor>(R() * Ainv));
    auto &D = *gibbs_D_;

    if (not gibbs_r_Rbeta_) gibbs_r_Rbeta_.reset(new VectorXd(r() - R() * beta_));
    auto &r_minus_Rbeta_ = *gibbs_r_Rbeta_;

    double s = std::sqrt(s2_);

    if (not gibbs_last_z_) {
        // If we don't have an initial value, draw an *untruncated* value and give it to
        // gibbsInitialize() to fix up.

        for (int trial = 1; ; trial++) {
            try { gibbsInitialize(Linear::draw(), 10*restrict_size_); break; }
            catch (constraint_failure&) {
                if (trial >= 10) throw;
                // Otherwise try again
            }
        }

        // gibbs_last_*_ will be set up now (or else we threw an exception)
    }
    auto &z = *gibbs_last_z_;
    auto &sigma = gibbs_last_sigma_;
    double s_sigma = 0;
    auto &rng = Random::rng();

    int num_draws = 1;
    if (gibbs_draws_ < draw_gibbs_burnin)
        num_draws += draw_gibbs_burnin - gibbs_draws_;
    else if (draw_gibbs_thinning > 1)
        num_draws = draw_gibbs_thinning;

    std::normal_distribution<double> stdnorm(0, 1);
    std::chi_squared_distribution<double> chisq(n_);
    boost::math::normal_distribution<double> stdnorm_dist(stdnorm.mean(), stdnorm.stddev());
    boost::math::chi_squared_distribution<double> chisq_dist(chisq.n());
    // Get this value now (if we're going to need it) to potentially avoid some extra cdf lookups in
    // the truncated distribution draw:
    const double chi_sq_median = restrict_size_ == 0 ? NAN : median(chisq_dist);

    for (int t = 0; t < num_draws; t++) { // num_draws > 1 if thinning or burning in

        // First take a sigma^2 draw that agrees with the previous z draw
        if (restrict_size_ == 0) {
            // If no restrictions, simple: just draw it from the unconditional distribution
            // NB: gamma(n/2,2) === chisq(v)
            sigma = std::sqrt(n_ / chisq(rng));
        }
        else {
            // Otherwise we need to look at the z values we drew above and draw and admissable
            // sigma^2 value from the range of values that wouldn't have caused a constraint
            // violation had we used it to form beta = beta_ + sigma*s*Ainv*z.  (See the method
            // documentation for the thorough details)
            double sigma_l = 0, sigma_u = INFINITY;
            for (size_t i = 0; i < restrict_size_; i++) {
                double denom = D.row(i) * z;
                if (denom == 0) {
                    // This means our z draw was exactly parallel to the restriction line--in which
                    // case, any amount of scaling will not violate the restriction (since we
                    // already know Z satisfies this restriction), so we don't need to do anything.
                    // (This case seems extremely unlikely, but just in case).
                    continue;
                }
                double limit = r_minus_Rbeta_[i] / (s * denom);

                if (denom > 0) {
                    if (limit < sigma_u) sigma_u = limit;
                }
                else {
                    if (limit > sigma_l) sigma_l = limit;
                }
            }

#ifdef ERIS_DEBUG
            try {
#endif

            // We want sigma s.t. beta_ + sigma * (s * A) * z has the right distribution for a
            // multivariate t, which is sigma ~ sqrt(n_ / chisq(n_)), *but* truncated to [sigma_l,
            // sigma_u].  To accomplish that truncation for sigma, we need to truncate the chisq(n_)
            // to [n_/sigma_u^2, n/sigma_l^2].
            sigma = std::sqrt(n_ / truncDist(chisq_dist, chisq, n_ / (sigma_u*sigma_u), n_ / (sigma_l*sigma_l), chi_sq_median));

#ifdef ERIS_DEBUG
            }
            catch (draw_failure &f) {
                throw draw_failure(f.what(), *this);
            }
#endif
        }

        s_sigma = sigma*s;

        for (unsigned int j = 0; j < K_; j++) {
            // Temporarily set the coefficient to 0, so that we don't have to maintain a bunch of
            // one-column-removed vectors and matrices below
            z[j] = 0.0;

            // Figure out l_j and u_j, the most binding constraints on z_j
            double lj = -INFINITY, uj = INFINITY;
            for (size_t r = 0; r < restrict_size_; r++) {
                // NB: not calculating the whole LHS vector and RHS vector at once, because there's
                // a good chance of 0's in the LHS vector, in which case we don't need to bother
                // calculting the RHS at all
                auto &dj = D(r, j);
                if (dj != 0) {
                    // Take the other z's as given, find the range for this one
                    double constraint = (r_minus_Rbeta_[r] / s_sigma - (D.row(r) * z)) / dj;
                    if (dj > 0) { // <= constraint (we didn't flip the sign by dividing by dj)
                        if (constraint < uj) uj = constraint;
                    }
                    else { // >= constraint
                        if (constraint > lj) lj = constraint;
                    }
                }
            }

            // Now lj is the most-binding bottom constraint, uj is the most-binding upper
            // constraint.  Make sure they aren't conflicting:
            if (lj >= uj) throw draw_failure("drawGibbs(): found impossible-to-satisfy linear constraints", *this);

#ifdef ERIS_DEBUG
            try {
#endif

            // Our new Z is a truncated standard normal (truncated by the bounds we just found)
            z[j] = truncDist(stdnorm_dist, stdnorm, lj, uj, 0.0);

#ifdef ERIS_DEBUG
            }
            catch (draw_failure &f) {
                throw draw_failure(f.what(), *this);
            }
#endif
        }

        gibbs_draws_++;
    }

    if (last_draw_.size() != K_ + 1) last_draw_.resize(K_ + 1);

    last_draw_.head(K_) = beta_ + s_sigma * Ainv * z;
    last_draw_[K_] = s_sigma * s_sigma;

    return last_draw_;
}

const VectorXd& LinearRestricted::drawRejection(long max_discards) {
    last_draw_mode = DrawMode::Gibbs;
    if (max_discards < 0) max_discards = draw_rejection_max_discards;
    draw_rejection_discards_last = 0;
    for (bool redraw = true; redraw; ) {
        redraw = false;
        auto &theta = Linear::draw();
        if (restrict_size_ > 0) {
            VectorXd Rbeta = restrict_select_.topRows(restrict_size_) * theta.head(K_);
            if ((Rbeta.array() > restrict_values_.head(restrict_size_).array()).any()) {
                // Restrictions violated
                redraw = true;
                ++draw_rejection_discards_last;
                ++draw_rejection_discards;
                if (draw_rejection_discards_last > max_discards) {
                    throw draw_failure("draw() failed: maximum number of inadmissible draws reached.");
                }
            }
        }
    }

    ++draw_rejection_success;

    return last_draw_;
}

double LinearRestricted::predict(const Eigen::Ref<const Eigen::RowVectorXd> &Xi) {
    return predict(Xi, 1000);
}

double LinearRestricted::predict(const Eigen::Ref<const Eigen::RowVectorXd> &Xi, long min_draws) {
    if (noninformative_)
        throw std::logic_error("Cannot call predict() on noninformative model");

    if (min_draws > mean_beta_draws_) {
        // First sum up new draws to make up the difference:
        VectorXd new_beta = VectorXd::Zero(K_);
        for (long i = mean_beta_draws_; i < min_draws; i++) {
            new_beta += draw().head(K_);
        }
        // Turn into a mean:
        new_beta /= min_draws;

        // If we had no draws at all before, just use the new vector
        if (mean_beta_draws_ == 0) {
            mean_beta_.swap(new_beta);
        }
        else {
            // Otherwise we need to combine means by calculating a weighted mean of the two means,
            // weighted by each mean's proportion of draws:
            double w_existing = (double) mean_beta_draws_ / min_draws;
            mean_beta_ = w_existing * mean_beta_ + (1 - w_existing) * new_beta;
        }

        mean_beta_draws_ = min_draws;
    }

    return Xi * mean_beta_;
}

bool LinearRestricted::hasRestriction(size_t k, bool upper) const {
    for (size_t row = 0; row < restrict_size_; row++) {
        // Only look at rows with a single non-zero coefficient:
        if ((restrict_select_.row(row).array() != 0).count() != 1)
            continue;

        const double &coef = restrict_select_(row, k);
        if (upper) {
            if (coef > 0) return true;
        }
        else {
            if (coef < 0) return true;
        }
    }
    return false;
}

double LinearRestricted::getRestriction(size_t k, bool upper) const {
    double most_binding = std::numeric_limits<double>::quiet_NaN();
    for (size_t row = 0; row < restrict_size_; row++) {
        // Only look at rows with a single non-zero coefficient:
        if (restrict_select_.row(row).cwiseEqual(0).count() != K()-1)
            continue;

        const double &coef = restrict_select_(row, k);
        if (coef == 0) continue;
        if (upper) {
            if (coef > 0) {
                double r = restrict_values_[row] / coef;
                if (std::isnan(most_binding) or r < most_binding) most_binding = r;
            }
        }
        else {
            if (coef < 0) {
                double r = restrict_values_[row] / coef;
                if (std::isnan(most_binding) or r > most_binding) most_binding = r;
            }
        }
    }
    return most_binding;
}

LinearRestricted::RestrictionProxy::RestrictionProxy(LinearRestricted &lr, size_t k, bool upper)
    : lr_(lr), k_(k), upper_(upper)
{}

bool LinearRestricted::RestrictionProxy::restricted() const {
    return lr_.hasRestriction(k_, upper_);
}
LinearRestricted::RestrictionProxy::operator double() const {
    return lr_.getRestriction(k_, upper_);
}

LinearRestricted::RestrictionProxy& LinearRestricted::RestrictionProxy::operator=(double r) {
    double Rk = upper_ ? 1.0 : -1.0;
    if (not upper_) r = -r;

    RowVectorXd R = RowVectorXd::Zero(lr_.K());
    R[k_] = Rk;

    lr_.addRestriction(R, r);

    return *this;
}

LinearRestricted::RestrictionIneqProxy::RestrictionIneqProxy(LinearRestricted &lr, size_t k)
    : lr_(lr), k_(k)
{}

bool LinearRestricted::RestrictionIneqProxy::hasUpperBound() const {
    return lr_.hasRestriction(k_, true);
}

double LinearRestricted::RestrictionIneqProxy::upperBound() const {
    return lr_.getRestriction(k_, true);
}

bool LinearRestricted::RestrictionIneqProxy::hasLowerBound() const {
    return lr_.hasRestriction(k_, false);
}

double LinearRestricted::RestrictionIneqProxy::lowerBound() const {
    return lr_.getRestriction(k_, false);
}

LinearRestricted::RestrictionIneqProxy& LinearRestricted::RestrictionIneqProxy::operator<=(double r) {
    RowVectorXd R = RowVectorXd::Zero(lr_.K());
    R[k_] = 1.0;
    lr_.addRestriction(R, r);
    return *this;
}

LinearRestricted::RestrictionIneqProxy& LinearRestricted::RestrictionIneqProxy::operator>=(double r) {
    RowVectorXd R = RowVectorXd::Zero(lr_.K());
    R[k_] = -1.0;
    lr_.addRestriction(R, -r);
    return *this;
}

LinearRestricted::operator std::string() const {
    std::ostringstream os;
    os << Linear::operator std::string();
    if (restrict_size_ == 0) os << "  No restrictions.\n";
    else if (restrict_size_ == 1) os << "  1 restriction:\n";
    else os << "  " << restrict_size_ << " restrictions:\n";
    for (size_t r = 0; r < restrict_size_; r++) {
        bool first = true;
        bool negate = (restrict_select_.row(r).array() <= 0).all();
        for (unsigned int j = 0; j < K_; j++) {
            double d = restrict_select_(r, j);
            if (negate) d = -d;
            if (d != 0) {
                if (first) {
                    first = false;
                    os << "    ";
                    if (d < 0) { os << "-"; d = -d; }
                }
                else {
                    if (d < 0) { os << " - "; d = -d; }
                    else       { os << " + "; }
                }
                if (d != 1) os << d << " ";
                os << "beta[" << j << "]";
            }
        }
        // If there were no non-zero elements then this is a trivial restriction of 0 <= r.
        if (first) { os << "    0"; negate = false; }

        os << (negate ? u8" ⩾ " : u8" ⩽ ") << (negate ? -1 : 1) * restrict_values_[r] << "\n";
    }

    return os.str();
}

std::string LinearRestricted::display_name() const { return "LinearRestricted"; }

}}
