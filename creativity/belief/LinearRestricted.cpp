// Eigen uses deprecated binder1st and binder2nd; until that is fixed, ignore the generated warning:
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include "creativity/belief/LinearRestricted.hpp"
#include <cmath>
#include <Eigen/QR>
#include <eris/algorithms.hpp>
#include <eris/Random.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/gamma.hpp>

#include <eris/debug.hpp>

using namespace Eigen;
using eris::Random;
using boost::math::erfc_inv;
using boost::math::gamma_p;
using boost::math::gamma_q;
using boost::math::gamma_p_inv;
using boost::math::gamma_q_inv;

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
    gibbs_alpha_.reset();
    gibbs_last_.reset();
    gibbs_draws_ = 0;
}

Ref<const MatrixXd> LinearRestricted::R() const {
    return restrict_select_.topRows(restrict_size_);
}

Ref<const VectorXd> LinearRestricted::r() const {
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
        if (not gibbs_last_ or gibbs_last_->size() != K_+1) gibbs_last_.reset(new VectorXd(K_+1));
        gibbs_last_->head(K_) = A * initial.head(K_);
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
                gibbs_last_.reset();
                throw constraint_failure("gibbsInitialize() couldn't find a way to satisfy the model constraints");
            }

            // Select a random constraint to fix:
            size_t fix = violations[std::uniform_int_distribution<size_t>(0, violations.size()-1)(rng)];

            // Aim at the nearest point on the boundary and overshoot (by 50%):
            adjusted += -overshoot * v[fix] / restrict_select_.row(fix).squaredNorm() * restrict_select_.row(fix).transpose();
        }

        if (not gibbs_last_ or gibbs_last_->size() != K_+1) gibbs_last_.reset(new VectorXd(K_+1));
        gibbs_last_->head(K_) = A * adjusted;
    }

    // Draw an appropriate sigma^2 value from the unconditional distribution, which should be
    // suitable as an initial value (subsequent draws will be from the conditional distribution
    // based on drawn beta values)
    (*gibbs_last_)[K_] = 1.0 / std::gamma_distribution<double>(n_/2, 2/(s2_*n_))(rng);
}

VectorXd LinearRestricted::gibbsLast() {
    return gibbs_last_ and gibbs_last_->size() == K_+1
        ? VcholL() * gibbs_last_->head(K_)
        : VectorXd();
}

const VectorXd& LinearRestricted::drawGibbs() {
    last_draw_mode = DrawMode::Gibbs;

    const MatrixXd &A = VcholLinv();
    const MatrixXd &Ainv = VcholL();

    if (not gibbs_D_) gibbs_D_.reset(
            restrict_size_ == 0
            ? new MatrixXd(0, K_) // It's rather pointless to use drawGibbs() with no restrictions, but allow it for debugging purposes
            : new MatrixXd(restrict_select_.topRows(restrict_size_) * Ainv));
    auto &D = *gibbs_D_;

    int num_draws = 1;
    if (gibbs_draws_ < draw_gibbs_burnin)
        num_draws = draw_gibbs_burnin - gibbs_draws_ + 1;
    else if (draw_gibbs_thinning > 1)
        num_draws = draw_gibbs_thinning;

    if (not gibbs_alpha_) gibbs_alpha_.reset(new VectorXd(A * beta_));
    auto &alpha = *gibbs_alpha_;

    auto &rng = Random::rng();
    if (not gibbs_last_) {
        // If we don't have an initial value, draw an *untruncated* value and give it to
        // gibbsInitialize() to fix up.
        
        for (int trial = 1; ; trial++) {
            try { gibbsInitialize(Linear::draw(), 10*restrict_size_); break; }
            catch (constraint_failure&) {
                if (trial >= 10) throw;
                // Otherwise try again
            }
        }

        // gibbs_last_ will be set up now (or else we threw an exception)
    }

    auto &Z = *gibbs_last_;


    for (int t = 0; t < num_draws; t++) { // num_draws > 1 if thinning or burning in

        for (unsigned int j = 0; j < K_; j++) {
            // Temporarily set the coefficient to 0, so that we don't have to maintain a bunch of
            // D-with-one-column-removed matrices below.
            Z[j] = 0.0;

            // Figure out l_j and u_j, the most binding constraints on \beta_j
            double lj = -INFINITY, uj = INFINITY;
            for (size_t r = 0; r < restrict_size_; r++) {
                // NB: not calculating the whole LHS vector and RHS vector at once, because there's
                // a good chance of 0's in the LHS vector, in which case we don't need to bother
                // calculting the RHS at all
                auto &dj = D(r, j);
                if (dj != 0) {
                    // NB: in the paper, the last term here is D_{-j} * z_{-j}, but because z_{j} =
                    // 0 is set up, this is equivalent (it involves multiplying and adding a zero,
                    // but that's less computational work than screwing around with removing
                    // elements from matrices, and much less programming work).
                    double constraint = (restrict_values_[r] - D.row(r) * Z.head(K_)) / dj;
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
            if (lj >= uj) throw draw_failure("drawGibbs(): found impossible-to-satisfy linear constraints");

            // Our new Z is a truncated normal (truncated by the bounds we just found)
            Z[j] = truncNorm(alpha[j], std::sqrt(Z[K_]), lj, uj);
        }

        // Now we need a sigma^2 draw.
        if (restrict_size_ == 0) {
            // If no restrictions, simple, just draw it from the unconditional distribution
            Z[K_] = 1.0 / std::gamma_distribution<double>(n_/2, 2/(s2_*n_))(rng);
        }
        else {
            // Otherwise we need to look at the Z values we drew above and draw and admissable
            // sigma^2 value from the range of values that wouldn't have caused a constraint
            // violation had we used it instead of Z[K_], then draw a new Z[K_] from that truncated
            // distribution.  (See the method documentation for the thorough details)
            double sl = 0, su = INFINITY;
            VectorXd Y = (Z.head(K_) - alpha) / std::sqrt(Z[K_]);
            for (size_t i = 0; i < restrict_size_; i++) {
                MatrixXd RiAinv = restrict_select_.row(i) * Ainv;
                double denom = (RiAinv * Y)(0,0);
                if (denom == 0) {
                    // This means our Y draw was exactly parallel to the restriction line--in which
                    // case, any amount of scaling will not violate the restriction (since we
                    // already know Z satisfies this restriction), so we don't need to do anything.
                    // This case seems extremely unlikely.
                    continue;
                }

                double numer = restrict_values_[i] - (RiAinv * alpha)(0,0);

                if (denom > 0)
                    su = std::min(su, numer / denom);
                else
                    sl = std::max(sl, numer / denom);
            }

            Z[K_] = 1.0 / truncGamma(n_/2, 2/(s2_*n_), 1.0/(su*su), 1.0/(sl*sl));
        }

        gibbs_draws_++;
    }

    if (last_draw_.size() != K_ + 1) last_draw_.resize(K_ + 1);

    last_draw_[K_] = Z[K_];
    last_draw_.head(K_) = Ainv * Z.head(K_);

    return last_draw_;
}

double LinearRestricted::truncDist(
                double min,
                double max,
                const std::function<double(double)> &cdf,
                const std::function<double(double)> &cdf_complement,
                const std::function<double(double)> &quantile,
                const std::function<double(double)> &quantile_complement,
                double median) {

    if (min >= max) throw std::runtime_error("Can't call truncDist() with min >= max!");

    double alpha, omega;
    bool alpha_comp, omega_comp;
    if (std::isnan(median)) { // No median, so use an extra call to figure out if we need complements
        alpha = cdf(min);
        alpha_comp = alpha > 0.5;
        if (alpha_comp) {
            alpha = cdf_complement(min);
            // min is right of median, so max will be too
            omega_comp = true;
            omega = cdf_complement(max);
        }
        else {
            // min was left of median.  Let's guess that max is right of the median first, which
            // seems like it might be slightly more likely and so might avoid a second cdf inverse
            // more often.
            omega = cdf_complement(max);
            omega_comp = omega < 0.5;
            if (not omega_comp) omega = cdf(max);
        }
    }
    else {
        alpha_comp = min > median;
        omega_comp = max > median;
        alpha = alpha_comp ? cdf_complement(min) : cdf(min);
        omega = omega_comp ? cdf_complement(max) : cdf(max);
    }

    if (not alpha_comp and omega_comp) {
        // min is left of the median and max is right: we need either make both complements or both
        // non-complements.  The most precision is going to be lost by whichever value is smaller,
        // so complement the larger value (so we'll either get both complements, or both
        // non-complements).
        if (alpha > omega) { alpha = 1 - alpha; alpha_comp = true; }
        else { omega = 1 - omega; omega_comp = false; }
    }
    // (alpha_comp and not omega_comp) is impossible, because that would mean 'min > max' and we
    // would have thrown and exception above.

    if (alpha_comp) {
        // Check for underflow (essentially: is 1-alpha equal to or closer to 0 than a double can
        // represent without reduced precision)
        if (alpha == 0 or std::fpclassify(alpha) == FP_SUBNORMAL)
            throw std::underflow_error("LinearRestricted::truncDist(): Unable to draw from truncated distribution: truncation range is too far in the upper tail");

        // Both alpha and omega are complements, so take a draw from [omega,alpha], then pass it
        // through the quantile_complement
        double u_comp = std::uniform_real_distribution<double>(omega, alpha)(Random::rng());
        return quantile_complement(u_comp);
    }
    else {
        // Check for underflow (essentially, is omega equal to or closer to 0 than a double can
        // represent without reduced precision)
        if (omega == 0 or std::fpclassify(omega) == FP_SUBNORMAL)
            throw std::underflow_error("LinearRestricted::truncDist(): Unable to draw from truncated distribution: truncation range is too far in the lower tail");

        // Otherwise they are ordinary cdf values, draw from the uniform and invert:
        double u = std::uniform_real_distribution<double>(alpha, omega)(Random::rng());
        return quantile(u);
    }
}

constexpr double sqrt2 = 1.41421356237309515;
double LinearRestricted::truncNorm(double mean, double sd, double min, double max) {
    if (min >= max) throw constraint_failure("Can't call truncNorm() with min >= max!");

    if (std::isnan(mean) or std::isnan(sd)) throw std::runtime_error("truncNorm called with NaN parameters");

    if (min < 0 and std::isinf(min) and max > 0 and std::isinf(max))
        return mean + sd * Random::rstdnorm();

    return truncDist(min, max,
            [&mean,&sd](double x) { return 0.5*std::erfc((mean-x) / (sd*sqrt2)); }, // cdf
            [&mean,&sd](double x) { return 0.5*std::erfc((x-mean) / (sd*sqrt2)); }, // cdf complement
            [&mean,&sd](double p) { return mean - sd*sqrt2 * erfc_inv(2*p); }, // cdf inverse
            [&mean,&sd](double q) { return mean + sd*sqrt2 * erfc_inv(2*q); }, // cdf complement inverse
            mean); // median (== mean)
}

double LinearRestricted::truncGamma(double shape, double scale, double min, double max) {
    if (min >= max) throw constraint_failure("Can't call truncGamma() with min >= max");
    if (max <= 0) throw constraint_failure("Can't call truncGamma() with non-positive max");

    if (min <= 0 and std::isinf(max))
        return std::gamma_distribution<double>(shape, scale)(Random::rng());

    return truncDist(min, max,
            // NB: gamma_p/gamma_q don't like being given infinity (but don't mind very large, but
            // finite, numbers); likewise the inverse ones give back DBL_MAX instead of INFINITY if
            // given 1 (or 0, for complement), so handle give both cases specially (the other end of
            // the distribution, where x=0, gives back 0/1 as expected).
            [&shape,&scale](double x) { return x > 0 and std::isinf(x) ? 1.0 : gamma_p(shape, x / scale); }, // cdf
            [&shape,&scale](double x) { return x > 0 and std::isinf(x) ? 0.0 : gamma_q(shape, x / scale); }, // cdf complement
            [&shape,&scale](double p) { return p == 1 ? INFINITY : scale * gamma_p_inv(shape, p); }, // cdf inverse
            [&shape,&scale](double q) { return q == 0 ? INFINITY : scale * gamma_q_inv(shape, q); }); // cdf complement inverse
}


/* This is the beginning of Geweke (1991).  Discarded.
    MatrixXd R_lindep(restrict_size_, K_);
    int R_lindep_used = 0;
    MatrixXd R = MatrixXd::Zero(K_, K_);
    VectorXd r_a = VectorXd::Constant(K_, -std::numeric_limits<double>::infinity());
    VectorXd r_b = VectorXd::Constant(K_, std::numeric_limits<double>::infinity());
    unsigned int last_rank = 0;
    ColPivHouseholderQR<MatrixXd> qr;
    for (unsigned int i = 0; i < restrict_size_ and last_rank < K_; i++) {
        // Try adding the row
        R.row(last_rank) = restrict_select_.row(i);
        qr = R.colPivHouseholderQr();
        unsigned int new_rank = qr.rank();
        if (new_rank == last_rank + 1) {
            // Great, adding the row increased the rank, keep it and copy the restricted value:
            r_b[last_rank] = restrict_values_[i];
            last_rank = new_rank;
        }
        else if (new_rank == last_rank) {
            // Otherwise the rank didn't increase, and there are two possibilities: this could be an
            // opposite bound of an earlier restriction, or it could be some more complicated linear
            // restriction.  The former we can handle (since the procedure allows two bounds for
            // each row); the latter we'll have to throw away and check after we draw.
            bool found_equiv = false;
            for (unsigned int prev = 0; prev < last_rank; prev++) {
                // NB: we only handle exact negatives; something like R_j = -2R_i won't be caught here and will just end up in the "can't handle restrictions" set.
                if (R.row(prev) == -restrict_select_.row(i) and r_a[prev] == -std::numeric_limits<double>::infinity()) {
                    r_a[prev] = -restrict_values_[i];
                    found_equiv = true;
                    break;
                }
            }
            if (!found_equiv) {
                R_lindep.row(R_lindep_used) = restrict_select_.row(i);
                R_lindep_used++;
            }

            // NB: R still has the problem row in it, but it'll be overwritten by the next pass (or
            // below)
        }
        else {
            // Safety check
            throw std::runtime_error("Internal LinearRestricted error: adding one row changed rank from " + std::to_string(last_rank) + " to " + std::to_string(new_rank) + ", which shouldn't be possible!");
        }
    }

    // If we don't have full rank, we'll have to fill out the rest of R by adding new rows of of all
    // zeros except for 1 column, and keeping any that increase the rank.  The associated values,
    // however, are left unrestricted because these aren't meant to be binding restrictions.
    unsigned int try_col = 0;
    while (last_rank < K_ and try_col < K_) {
        // Add a new row.
        for (unsigned int i = 0; i < K_; i++) {
            R(last_rank, i) = (i == try_col ? 1 : 0);
        }

        qr = R.colPivHouseholderQr();
        unsigned int new_rank = qr.rank();
        if (new_rank == last_rank + 1) {
            last_rank = new_rank;
        }
        else if (new_rank != last_rank) {
            // Safety check
            throw std::runtime_error("Internal LinearRestricted error: adding one row changed rank from " + std::to_string(last_rank) + " to " + std::to_string(new_rank) + ", which shouldn't be possible!");
        }
    }

    if (last_rank < K_) {
        // Don't throw a draw_failure here, this is more serious and should be impossible.
        throw std::runtime_error("Internal LinearRestricted error: failed to make R full-rank, unable to draw()");
    }

    // Now we've got an invertible R, hurray!
    MatrixXd R_inv = qr.inverse();
    */

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

void LinearRestricted::print(std::ostream &os) const {
    Linear::print(os);
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
}

std::string LinearRestricted::print_name() const { return "LinearRestricted"; }

}}
