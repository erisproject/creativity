// Eigen uses deprecated binder1st and binder2nd; until that is fixed, ignore the generated warning:
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include "creativity/belief/LinearRestricted.hpp"
#include <cmath>

#include <eris/debug.hpp>

using namespace Eigen;

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
}

void LinearRestricted::addRestrictionsGE(const Ref<const MatrixXd> &R, const Ref<const VectorXd> &r) {
    return addRestrictions(-R, -r);
}

void LinearRestricted::clearRestrictions() {
    restrict_size_ = 0;
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
    draw_discards = 0;
    draw_discards_cumulative = 0;
    draw_success_cumulative = 0;
}

const Ref<const MatrixXd> LinearRestricted::R() const {
    return restrict_select_.topRows(restrict_size_);
}

const Ref<const VectorXd> LinearRestricted::r() const {
    return restrict_values_.head(restrict_size_);
}

const VectorXd& LinearRestricted::draw() {
    bool redraw = true;
    draw_discards = 0;
    while (redraw) {
        redraw = false;
        auto &theta = Linear::draw();
        /*
        ERIS_DBGVAR(restrict_size_);
        if (restrict_size_ > 0) {
            ERIS_DBGVAR(theta.head(K_).transpose());
            ERIS_DBGVAR(restrict_select_);
            ERIS_DBGVAR((restrict_select_.topRows(restrict_size_) * theta.head(K_)).array());
            ERIS_DBGVAR(restrict_values_.head(restrict_size_));
        }
        */
        if (restrict_size_ > 0 and not (
                (
                    (restrict_select_.topRows(restrict_size_) * theta.head(K_)).array()
                        // Rbeta
                    <=  // <=
                        // r
                    restrict_values_.head(restrict_size_).array()
                ).all()
        )) {
            redraw = true;
            ++draw_discards;
            ++draw_discards_cumulative;
            if (draw_discards > draw_discards_max) {
                throw draw_failure("draw() failed: maximum number of inadmissible draws reached.");
            }
        }
    }

    ++draw_success_cumulative;

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


}}
