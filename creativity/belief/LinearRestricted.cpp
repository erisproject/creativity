#include "creativity/belief/LinearRestricted.hpp"
#include <cmath>

using namespace Eigen;

namespace creativity { namespace belief {

std::vector<double>::iterator LinearRestricted::lowerBounds() {
    if (restrict_ge_.size() <= K_) restrict_ge_.resize(K_+1, std::numeric_limits<double>::quiet_NaN());
    return restrict_ge_.begin();
}

std::vector<double>::iterator LinearRestricted::upperBounds() {
    if (restrict_le_.size() <= K_) restrict_le_.resize(K_+1, std::numeric_limits<double>::quiet_NaN());
    return restrict_le_.begin();
}

void LinearRestricted::addRestriction(RowVectorXd R, double r) {
    restrict_linear_.emplace_back(std::move(R), std::move(r));
}

void LinearRestricted::clearRestrictions() {
    restrict_ge_.clear();
    restrict_le_.clear();
    restrict_linear_.clear();
}

void LinearRestricted::discard(unsigned int burn) {
    Linear::discard(burn); // Pass it up (mainly for the no-empty-model check)
    mean_beta_draws_ = 0; // Reset the number of beta draws so that the next predict() redraws
}

void LinearRestricted::discardForce(unsigned int burn) {
    Linear::discardForce(burn);
    mean_beta_draws_ = 0;
}

const VectorXd& LinearRestricted::draw() {
    bool redraw = true;
    draw_discards = 0;
    while (redraw) {
        redraw = false;
        auto &theta = Linear::draw();
        // First check simple, one-parameter <=/>= restrictions:
        if (not restrict_le_.empty() or not restrict_ge_.empty()) {
            for (unsigned int i = 0; i < K_; i++) {
                if (
                        (not restrict_le_.empty() and std::isfinite(restrict_le_[i]) and not(theta[i] <= restrict_le_[i]))
                        or
                        (not restrict_ge_.empty() and std::isfinite(restrict_ge_[i]) and not(theta[i] >= restrict_ge_[i]))
                ) {
                    // beta[i] violates one of its restrictions; try again
                    redraw = true;
                    break;
                }
            }
        }
        // If no single-parameter restriction was violated, now check any more general linear
        // restrictions:
        if (not redraw and not restrict_linear_.empty()) {
            for (auto &r : restrict_linear_) {
                if (not (r.first * theta.head(K_) <= r.second)) {
                    redraw = true;
                    break;
                }
            }
        }

        if (redraw) {
            ++draw_discards;
            ++draw_discards_cumulative;
            if (draw_discards > draw_discards_max) {
                throw std::runtime_error("draw() failed: maximum number of inadmissible draws reached.");
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


}}
