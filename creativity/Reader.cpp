#include "creativity/Creativity.hpp"
#include "creativity/Reader.hpp"
#include "creativity/Book.hpp"
#include "creativity/BookMarket.hpp"
#include "creativity/CopyrightPolice.hpp"
#include "creativity/belief/Profit.hpp"
#include <eris/debug.hpp>
#include <eris/random/util.hpp>
#include <boost/random/chi_squared_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <Eigen/Core>
#include <algorithm>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <map>
#include <list>
#include <cmath>

using namespace eris;
using namespace Eigen;
using namespace eris::learning;
using namespace creativity::belief;

namespace creativity {

constexpr std::initializer_list<double> Reader::default_distance_penalty_polynomial;
constexpr std::initializer_list<double> Reader::default_num_books_penalty_polynomial;
const std::vector<unsigned int> Reader::profit_stream_ages{{1,2,4,8}};

Reader::Reader(Creativity &creativity, const Position &pos) :
    WrappedPositional<Agent>(pos, creativity.parameters.boundary, -creativity.parameters.boundary),
    creativity_{creativity},
    profit_belief_{new Profit()},
    profit_belief_extrap_{profit_belief_}
{
    profit_stream_beliefs_.emplace(1, ProfitStream(1));
    distancePenaltyPolynomial(default_distance_penalty_polynomial);
    numBooksPenaltyPolynomial(default_num_books_penalty_polynomial);
}

void Reader::distancePenaltyPolynomial(std::vector<double> coef) {
    // Trim off any trailing 0's to simplify the polynomial expression:
    while (not coef.empty() and coef.back() == 0) coef.pop_back();

    if (coef.size() > 1 and coef.back() < 0) throw std::domain_error("Invalid distancePenaltyPolynomial: last non-zero, non-constant coefficient must be positive.");

    // Make sure all coefficients are finite
    for (size_t i = 0; i < coef.size(); i++) {
        if (not std::isfinite(coef[i]))
            throw std::domain_error("Invalid uPolynomial: coef[" + std::to_string(i) + "] is not finite (" + std::to_string(coef[i]) + ")");
    }

    // First non-zero, non-constant coefficient must be positive
    for (size_t i = 1; i < coef.size(); i++) {
        if (coef[i] < 0) throw std::domain_error("Invalid distancePenaltyPolynomial: first non-zero, non-constant coefficient must be positive.");
        else if (coef[i] > 0) break;
    }

    double xlast = 0;
    double plast = Creativity::evalPolynomial(0, coef);
    for (const double &x : {1e-100, 1e-10, .1, 1., 10., 100., 1e10}) {
        double p = Creativity::evalPolynomial(x, coef);
        if (p < plast)
            throw std::domain_error("Invalid uPolynomial: polynomial is not increasing: f(" + std::to_string(x) + ") = " + std::to_string(p) +
                    " is less than f(" + std::to_string(xlast) + ") = " + std::to_string(plast));
        xlast = x;
        plast = p;
    }

    dist_penalty_poly_ = std::move(coef);
}

const std::vector<double>& Reader::distancePenaltyPolynomial() const {
    return dist_penalty_poly_;
}

double Reader::distancePenalty(double distance) const {
    return Creativity::evalPolynomial(distance, distancePenaltyPolynomial());
}

double Reader::u(double money, const std::unordered_map<eris::SharedMember<Book>, std::reference_wrapper<BookCopy>> &books) const {
    double u = money - numBooksPenalty(books.size());
    for (auto &bc : books) {
        u += uBook(bc.first, bc.second.get().quality);
    }
    return u;
}
const double& Reader::u() const {
    return u_curr_;
}
const double& Reader::uLifetime() const {
    return u_lifetime_;
}
const std::unordered_map<SharedMember<Book>, std::reference_wrapper<BookCopy>>& Reader::newBooks() const {
    return library_new_;
}
const std::set<SharedMember<Book>>& Reader::wrote() const {
    return wrote_;
}

void Reader::numBooksPenaltyPolynomial(std::vector<double> coef) {
    unsigned xlast = 0;
    double plast = Creativity::evalPolynomial(xlast, coef);
    for (const unsigned &x : {1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,1000,1000000}) {
        double p = Creativity::evalPolynomial(x, coef);
        if (p < plast)
            throw std::domain_error("Invalid numBooksPenaltyPolynomial: f(" + std::to_string(x) + ") = " + std::to_string(p) +
                    " is less than f(" + std::to_string(xlast) + ") = " + std::to_string(plast));
        plast = p;
        xlast = x;
    }

    nbooks_penalty_poly_ = std::move(coef);
}
const std::vector<double>& Reader::numBooksPenaltyPolynomial() const {
    return nbooks_penalty_poly_;
}
double Reader::numBooksPenalty(unsigned long n) const {
    return Creativity::evalPolynomial(n, numBooksPenaltyPolynomial());
}

double Reader::creationQuality(double effort) const {
    return creationQuality(creation_shape, creation_scale, effort);
}

double Reader::creationQuality(double shape, double scale, double effort) {
    if (effort < 0)
        throw std::domain_error("Reader::creationQuality() error: effort cannot be negative");
    if (shape >= 1)
        throw std::logic_error("Reader::creationQuality() error: shape must be less than 1");
    if (scale < 0)
        throw std::logic_error("Reader::creationQuality() error: scale cannot be negative");

    if (scale == 0) return 0.0; // This reader always produces 0 quality

    if (effort == 0) return 0.0; // All of the functional forms of the effort yield q(0) = 0

    return scale * (
            shape == 0.0 ? std::log(effort + 1.0) : // Special case for log
            shape == 0.5 ? 2.0 * (std::sqrt(effort + 1.0) - 1.0) : // Slightly more efficient sqrt calculation
            (std::pow(effort + 1.0, shape) - 1.0) / shape // Otherwise use the full formula
    );
}

double Reader::creationEffort(double quality) const {
    return creationEffort(creation_shape, creation_scale, quality);
}

double Reader::creationEffort(double shape, double scale, double quality) {
    if (quality < 0)
        throw std::domain_error("Reader::creationEffort() error: quality cannot be negative");
    if (shape >= 1)
        throw std::logic_error("Reader::creationEffort() error: shape cannot be negative");
    if (scale < 0)
        throw std::logic_error("Reader::creationEffort() error: scale cannot be negative");

    if (quality == 0) return 0.0; // All functional forms yield q(0) = 0, and so q^-1(0) = 0.
    else if (scale == 0) {
        // If scale is 0, but quality is not, the function is not invertible (because any
        // level of effort yields 0 quality), so return infinity.
        return std::numeric_limits<double>::infinity();
    }
    else if (scale < 0 and quality >= (scale / -shape)) {
        // When beta (scale) is below 0, there's an asymptote at alpha/-beta; no effort
        // level exists that can produce higher quality than that asymptote, so return infinity.
        return std::numeric_limits<double>::infinity();
    }

    return
        shape == 0.0 ? std::exp(quality / scale) - 1.0 : // The special log case
        std::pow(quality * shape / scale + 1.0, 1.0/shape) - 1.0; // Non-log case
}

double Reader::piracyCost() const {
    return creativity_.parameters.cost_piracy;
}

const std::unordered_set<SharedMember<Reader>>& Reader::friends() const {
    return friends_;
}

bool Reader::addFriend(SharedMember<Reader> new_pal, bool recurse) {
    auto inserted = friends_.insert(std::move(new_pal));
    if (inserted.second and recurse) (*inserted.first)->addFriend(sharedSelf(), false);
    return inserted.second;
}

bool Reader::removeFriend(const SharedMember<Reader> &old_pal, bool recurse) {
    auto found = friends_.find(old_pal);
    if (found == friends_.end()) return false; // Not found

    if (recurse) (*found)->removeFriend(sharedSelf(), false);
    friends_.erase(found);
    return true;
}


void Reader::interBegin() {
    // Move a random distance in a random direction
    if (creativity_.parameters.reader_step_mean > 0) {
        double step_dist = boost::random::chi_squared_distribution<double>()(random::rng()) * creativity_.parameters.reader_step_mean;
        if (step_dist != 0)
            moveBy(step_dist * Position::random(position().dimensions));
    }
}

void Reader::interOptimize() {
    // Update the various beliefs
    updateBeliefs();

    double income_available = assets[creativity_.money] + creativity_.disposableIncome();

    auto sim = simulation();
    const double &market_books_avg = creativity_.market_books_avg;
    const double authored_books = wrote().size();

    // Any books not in new_prices_ but that were on the market last period will be removed from the
    // market
    new_prices_.clear();

    auto &rng = random::rng();
    const double &cost_market = creativity_.parameters.cost_market,
          &cost_unit = creativity_.parameters.cost_unit;

    if (not wrote_market_.empty()) {
        if (usableBelief(demand_belief_)) {
            // Figure out the expected profitability of each book; positive ones will be kept on the market
            std::map<double, std::vector<std::pair<SharedMember<Book>, double>>, std::greater<double>> profitability;
            for (auto &book : wrote_market_) {
                auto book_mkt = book->market();

                try {
                    auto max = demand_belief_.argmaxP(creativity_.parameters.prediction_draws, book->qualityMean(), book->lifeSales(), book->lifePirated(), creativity_.parameters.readers,
                            sim->t() - 1 - book->lastSale(), book->age(), authored_books - 1, market_books_avg, cost_unit,
                            creativity_.parameters.income / 10);
                    const double &p = max.first;
                    const double &q = max.second;

                    const double profit = (p - cost_unit) * q - cost_market;

                    if (profit > 0) {
                        // Profitable to keep this book on the market, so do so
                        profitability[profit].push_back(std::make_pair(book, p));
                    }
                }
                catch (BayesianLinear::draw_failure &fail) {
                    // If we fail to get an admissable draw, pull the book from the market.
                }
            }

            // Store the decided prices for each book staying on the market in decreasing order of
            // profitability UNLESS we can't afford the fixed cost of keeping the book on the market.
            // Anything with a negative expected profit or costs that can't be covered will be removed from
            // the market.  Note that the costs calculated here aren't actually incurred until interApply().
            for (auto prof_iter = profitability.begin(), prof_iter_end = profitability.end();
                    income_available >= cost_market and prof_iter != prof_iter_end;
                    prof_iter++) {
                auto &books = prof_iter->second;
                // If we don't have enough to put all the books at this profitability level on the
                // market, shuffle them so that we choose randomly.
                if (income_available < cost_market * books.size())
                    std::shuffle(books.begin(), books.end(), rng);

                for (auto &b : books) {
                    income_available -= cost_market;
                    new_prices_.emplace(b.first, b.second);
                    // Stop if we run out of income:
                    if (income_available < cost_market) break;
                }
            }
        }
        // Otherwise bypass market predictions if our market demand belief isn't usable: use initial
        // model action parameters instead (unless we don't have enough income to support even a
        // single book, in which case we can't keep anything on the market).
        else if (creativity_.parameters.initial.prob_keep > 0 and income_available >= cost_market) {

            std::vector<SharedMember<Book>> books(wrote_market_.begin(), wrote_market_.end());
            if (income_available < books.size() * cost_market) {
                // If we don't have enough income to keep every book on the market (which is the
                // worst case scenario for income), shuffle the list so that we choose which ones to
                // keep randomly.
                std::shuffle(books.begin(), books.end(), rng);
            }
            for (auto &book : books) {
                if (random::rcoin(creativity_.parameters.initial.prob_keep)) {
                    income_available -= cost_market;
                    double new_price = (book->price() - cost_unit) * creativity_.parameters.initial.keep_price + cost_unit;
                    new_prices_.emplace(book, new_price);
                    // If we don't have enough income for any more books, leave the rest off the
                    // market:
                    if (income_available < cost_market) break;
                }
            }
        }
    }

    // If we have no book in the pipeline, we can consider whether or not to write (or write
    // randomly, if profit and demand beliefs aren't yet usable).
    if (create_countdown_ == -1) {

        // Make sure we have enough income to create a book, which requires (at 0 effort)
        // creation_fixed; if creation is instantaneous, we also need to have cost_market available
        // so that we can actually bring the book to market.
        double creation_base_cost = creativity_.parameters.creation_fixed;
        // Note: creation_time == 1 can (1/3 of the time) also mean the book is completed in the
        // current period, so make sure that we reserve enough income in case that happens
        if (creativity_.parameters.creation_time <= 1) creation_base_cost += cost_market;

        // Figure out what the minimum and maximum effort we can expend are.
        //
        // Set minimum effort to the effort required to get a book with a quality large enough to
        // cover a book at distance 0, which is a reader's first book, and is sold for a price 5%
        // above marginal cost.  The 5% is fairly innocuous: it's essentially just an extra amount
        // to reflect that distance will never be 0, and price will never fall all the way to
        // marginal cost.
        //
        // Marginal cost is the unit cost, but if we have public sharing, we'll use the public
        // sharing cost (which is the minimum of unit cost and piracy cost) if that is lower (since
        // that's what the public sharing marginal cost and price ends up being).
        double marginal_cost = creativity_.parameters.cost_unit;
        if (creativity_.publicSharingActive() and creativity_.parameters.cost_piracy < marginal_cost)
            marginal_cost = creativity_.parameters.cost_piracy;
        double min_effort = creationEffort(distancePenalty(0) + numBooksPenalty(1) + 1.05 * marginal_cost);
        double max_effort = income_available - creation_base_cost;
        // Don't even bother trying unless min_effort is actually obtainable
        if (income_available >= creation_base_cost and max_effort > min_effort) {
            // Create a book if maximized (wrt effort) `E(profit(effort)) - effort > 0`


            if (usableBelief(*profit_belief_) and usableBelief(demand_belief_)) {
                // NB: profit_belief_extrap_ is a copy or update of profit_belief_, so will also be usable

                // Find the l that maximizes creation profit given current wealth (which would have come
                // from past sales, if non-zero) plus the fixed income we're about to receive, minus
                // whatever we decided to spend above to keep books on the market.

                std::pair<double, double> effort_profit{-1.0, -1.0};
                create_price_ = std::numeric_limits<double>::quiet_NaN();

                try {
                    effort_profit = profit_belief_extrap_->argmaxL(
                            creativity_.parameters.prediction_draws,
                            [this] (double l) { return creationQuality(l); },
                            authored_books, market_books_avg,
                            min_effort, // l_min
                            income_available - creation_base_cost // l_max
                            );

                    // Also pick a price (just to make sure we *can* pick a price without draw
                    // failures).  If we create this period, we use it; if we don't we'll pick
                    // another one when we create (and fall back on this one if we can't draw then).

                    // As a safety against absurd prices, we also enforce a price ceiling of 20% of
                    // per-period income.
                    std::pair<double, double> max;
                    max = demand_belief_.argmaxP(creativity_.parameters.prediction_draws, create_quality_, 0, 0,
                            creativity_.parameters.readers, 0, 0, authored_books, market_books_avg,
                            cost_unit, creativity_.parameters.income / 10);
                    // If this is zero, that's okay: that might be a release-directly-to-public book
                    create_price_ = max.first;
                } catch (const BayesianLinear::draw_failure &e) {
                    ERIS_DBG("draw failure in profit_extrap argmaxL for reader=" << id() << ", t=" << simulation()->t() << ": " << e.what());
                    // If we get draw failures, ignore them (and don't create)
                }

                // If the optimal effort level gives a book that earns back at least the fixed
                // creation cost, we'll do it.
                if (effort_profit.first >= 0 and create_price_ > 0 and effort_profit.second > creation_base_cost) {
                    double quality = creationQuality(effort_profit.first);
                    create_quality_ = quality;
                    create_effort_ = effort_profit.first;
                    create_starting_ = true;
                    create_countdown_ = creativity_.parameters.creation_time;
                    create_position_ = position();
                }
            }
            else {
                // If we have no useful profit belief yet, just use the initial values:
                if (creativity_.parameters.initial.prob_write > 0 and random::rcoin(creativity_.parameters.initial.prob_write)) {
                    double effort = random::runiform(creativity_.parameters.initial.l_min,
                            creativity_.parameters.initial.l_min + creativity_.parameters.initial.l_range);
                    // Make sure the required effort doesn't exceed the available funds; if it does,
                    // reduce the effort:
                    if (creation_base_cost + effort > income_available)
                        effort = income_available - creation_base_cost;

                    create_starting_ = true;
                    create_countdown_ = creativity_.parameters.creation_time;
                    create_effort_ = effort;
                    create_quality_ = creationQuality(effort);
                    create_position_ = position();
                    create_price_ = creativity_.parameters.cost_unit +
                        boost::random::uniform_real_distribution<double>(
                                creativity_.parameters.initial.p_min,
                                creativity_.parameters.initial.p_min + creativity_.parameters.initial.p_range)(rng);
                }
            }

            if (create_countdown_ >= 1) {
                create_countdown_ += boost::random::uniform_int_distribution<int>(-1, 1)(random::rng());
            }
        }
    }

    if (create_countdown_ == 0 and !create_starting_) {
        // A currently-in-progress book is done and we didn't create it this period, so choose its
        // price based on current demand beliefs (which may have changed since the creation was
        // started).  If the optimal price is 0, we won't release it at all: the creation is a sunk
        // cost: it apparently doesn't seem profitable any more.  If we can't draw a price at all,
        // we'll just fall back on the initial one decided when the book was created.
        if (usableBelief(demand_belief_)) {
            // As a safety against absurd prices, we also enforce a price ceiling of 20% of
            // per-period income.
            std::pair<double, double> max;
            try {
                max = demand_belief_.argmaxP(creativity_.parameters.prediction_draws, create_quality_, 0, 0,
                        creativity_.parameters.readers, 0, 0, authored_books, market_books_avg,
                        cost_unit, creativity_.parameters.income / 10);
                create_price_ = max.first;
            } catch (BayesianLinear::draw_failure &e) {
                // Draw failure: don't change create_price_; we'll fall back to what it got set
                // to initially at creation time.
                ERIS_DBG("draw failure in demand argmaxP calculation; reader="<< id() << ", t=" << simulation()->t());
                ERIS_DBGVAR(demand_belief_.draw_rejection_success);
            }
        }
        else {
            // No usable beliefs; use initial behaviour parameters draw:
            create_price_ = creativity_.parameters.cost_unit +
                boost::random::uniform_real_distribution<double>(
                        creativity_.parameters.initial.p_min,
                        creativity_.parameters.initial.p_min + creativity_.parameters.initial.p_range)(rng);
        }
    }
}

void Reader::interApply() {
    // Give potential income (this has to be before authorship decision, since authorship requires
    // giving up some potential income: actual income will be this amount minus whatever is given up
    // to author)
    // NB: we give total income; any policy tax is removed by another policy agent, e.g.
    // PublicTracker, in a later-priority interApply().
    assets[creativity_.money] += creativity_.parameters.income;

    const Bundle cost_market(creativity_.money, creativity_.parameters.cost_market);

    auto sim = simulation();
    if (create_starting_) {
        // Incur the creation effort cost as soon as we start creating
        assets.transferApprox({ creativity_.money, create_effort_ + creativity_.parameters.creation_fixed }, 1e-6);
        create_starting_ = false;
    }
    if (create_countdown_ > 0) {
        create_countdown_--;
    }
    else if (create_countdown_ == 0) {
        // Book is finished, release it.

        if (create_price_ > 0) {
            // Author-selected price is positive, so let's release.

            // But first make sure we can afford it (keeping in mind that we'll have to pay our
            // taxes as well); if we can't, we'll just leave it here and (hopefully) release it next
            // period.
            if (assets.hasApprox(cost_market)) {
                // Remove the first period fixed cost of bringing the book to market:
                assets.transferApprox(cost_market, 1e-6);

                auto newbook = sim->spawn<Book>(creativity_, create_position_, sharedSelf(), wrote_.size(), create_quality_);
                sim->spawn<BookMarket>(creativity_, newbook, create_price_);

                /// If enabled, add some noise in a random direction to the position
                if (creativity_.parameters.book_distance_mean > 0) {
                    double step_dist = creativity_.parameters.book_distance_mean * boost::random::chi_squared_distribution<double>()(random::rng());
                    newbook->moveBy(step_dist * Position::random(create_position_.dimensions));
                }

                // NB: newbook can't go in a container yet; it will notify us (via register*
                // methods) when the book is added to the simulation and/or market.

                create_countdown_--;
            }
            // Otherwise we can't afford to bring it to market just now, so hold onto it until next
            // period (leave create_countdown_ at 0).
        }
        else {
            // According to demand, it isn't profitable to bring book to market, so just create it
            // but don't market it.  If we're in a pure-private, no-piracy environment, no one will
            // ever get the book; in a public provisioning world, there may be public payoff; in a
            // piracy world readers connected to the author can pirate it from the author.
            sim->spawn<Book>(creativity_, create_position_, sharedSelf(), wrote_.size(), create_quality_);
            // No market! (PublicTracker, if created, will figure this out and handle it).
            create_countdown_--;
        }
    }

    std::vector<SharedMember<Book>> remove;
    for (auto &b : wrote_market_) {
        // If it's staying on the market, update the price and incur the fixed cost
        if (new_prices_.count(b) > 0) {
            b->market()->setPrice(new_prices_[b]);
            // NB: we've already checked that we have enough income (if we didn't, we would not have
            // put it in new_prices_).
            assets.transferApprox(cost_market, 1e-6);
        }
        else {
            // No new price for an old book, which means we're removing the book from the market
            remove.push_back(b);
        }
    }

    for (auto &b : remove) {
        wrote_market_.erase(b);
        simulation()->remove(b->market());
    }
}

// Called when the book gets actually added to the simulation
void Reader::registerAuthoredBook(SharedMember<Book> book) {
    wrote_.insert(wrote_.end(), book);
    library_.emplace(SharedMember<Book>(book), BookCopy(create_quality_, BookCopy::Status::wrote, simulation()->t()));
}

// Calls when a market gets added or removed from the book
void Reader::registerMarketUpdate(SharedMember<Book> book) {
    if (book->hasPrivateMarket())
        wrote_market_.insert(book);
    else
        wrote_market_.erase(book);
}

double Reader::uBook(const SharedMember<Book> &b) const {
    return uBook(b, quality(b));
}

double Reader::uBook(const SharedMember<Book> &b, double quality) const {
    return quality - distancePenalty(distance(b));
}

double Reader::quality(const SharedMember<Book> &b) const {
    // If we already have a copy of this book, return the realized quality value
    auto found = library_.find(b);
    if (found != library_.end())
        return found->second.quality;

    // Otherwise return the mean quality
    return b->qualityMean();
}

const std::unordered_map<SharedMember<Book>, BookCopy>& Reader::library() const { return library_; }

const belief::Profit& Reader::profitBelief() const { return *profit_belief_; }
const belief::Profit& Reader::profitExtrapBelief() const { return *profit_belief_extrap_; }
bool Reader::profitExtrapBeliefDiffers() const { return profit_belief_ != profit_belief_extrap_; }
const belief::Demand& Reader::demandBelief() const { return demand_belief_; }
const belief::ProfitStream& Reader::profitStreamBelief(const unsigned int age, const bool usable) const {
    // Get the first belief > age
    auto it = profit_stream_beliefs_.upper_bound(age);
    if (it == profit_stream_beliefs_.begin())
        // This probably means age=0 got passed in, since profit_stream_beliefs_ starts with a 1.
        throw std::runtime_error("Invalid age (" + std::to_string(age) + ") passed to Reader::profitStreamBelief");
    // upper_bound gives first greater than the given value, so we need to back up one; and we may
    // need to back up more if only usable beliefs were requested
    do {
        it--;
    } while (usable and it != profit_stream_beliefs_.begin() and not usableBelief(it->second));

    return it->second;
}
const std::map<unsigned int, belief::ProfitStream>& Reader::profitStreamBeliefs() const {
    return profit_stream_beliefs_;
}

bool Reader::usableBelief(const BayesianLinear &model) const {
    return not model.noninformative() and model.n() - model.K() >= creativity_.parameters.initial.belief_threshold;
}

void Reader::updateBeliefs() {
    updateDemandBelief();
    //updateProfitStreamBelief();
    updateProfitBelief();

    // Clear any "on-market" book references that aren't on the market anymore, since they just got
    // incorporated into the beliefs above.
    std::list<SharedMember<Book>> remove;
    for (auto &book : book_cache_market_) {
        if (not book->hasPrivateMarket())
            remove.push_back(book);
    }
    for (auto &book : remove)
        book_cache_market_.erase(book);
}

void Reader::updateDemandBelief() {
    // If we aren't starting simulation period t=3 or later, don't do anything: we can't
    // incorporated books into the belief because we need the lagged market book count
    if (simulation()->t() < 3) return;

    double weaken = creativity_.priorWeight();

    // Figure out which books might yield usable data: i.e. those on the market
    std::list<SharedMember<Book>> mktbooks;
    for (auto &b : book_cache_market_) {
        if (b->hasPrivateMarket()) mktbooks.push_back(b);
    }

    if (not mktbooks.empty()) {
        VectorXd y(mktbooks.size());
        MatrixXd X(mktbooks.size(), demand_belief_.K());
        // NB: this runs in the interoptimizer, which means t has already been incremented
        auto last_t = simulation()->t() - 1;
        size_t i = 0;
        for (const auto &b : mktbooks) {
            y[i] = b->sales(last_t);
            X.row(i) = Demand::bookRow(b, b->qualityMean(), creativity_.market_books_avg);
            i++;
        }

        demand_belief_ = decltype(demand_belief_)(std::move(demand_belief_), y.head(i), X.topRows(i), weaken);
    }
    else if (weaken != 1.0)
        demand_belief_ = decltype(demand_belief_)(std::move(demand_belief_), weaken);

}
void Reader::updateProfitStreamBelief() {
    // FIXME: before re-enabling this, the belief needs to be updated somehow to track public prize
    // money
    throw std::runtime_error("Profit stream beliefs are not available");

    /*
    // Figure out which books might yield usable data
    std::list<SharedMember<Book>> potential;
    for (auto &b : book_cache_market_) {
        if (not b->hasMarket()) potential.push_back(b);
    }

    // Map book age into lists of books that survived on the market at least that long.  Only books
    // that have just left the market are considered (whether or not this reader bought or obtained
    // by piracy).
    std::map<unsigned int, std::vector<SharedMember<Book>>> to_learn;
    for (auto &book : potential) {
        if (not book->hasMarket()) { // The last time we checked, it had a market, so it just left
            unsigned long periods = book->marketPeriods();
            if (periods <= 1) continue; // We can't do anything with a single-period book

            // We only look for the predefined age values, and only include books with
            // marketPeriods() strictly greater than the predefined age values (because a book on
            // the market for x periods can only contribute to models with (x-1) past profit
            // variables.
            //
            // However, it's possible that the book was kept on the market for too long (i.e. it had
            // negative profits in the last 2+ periods), so reduce `periods` until either the last or
            // second-last referenced period had positive profits.  (The last having negative
            // profits is fine, because it will be the dependent variable for the model which only
            // goes up to the second-last, and we want the model to be able to predict the *first*
            // period the remaining profit becomes negative).
            //
            // Note that "profits" here are calculated as what profits would be with *this* reader's
            // fixed and unit costs, even if the actual author has different costs.  This is because
            // this belief is for predicting what the profit stream would have been for *this* reader.
            const unsigned long &created = book->created();
            double last_profit = book->revenue(created+periods-1) - cost_market - cost_unit*book->sales(created+periods-1);
            double second_last_profit = book->revenue(created+periods-2) - cost_market - cost_unit*book->sales(created+periods-2);
            while (periods > 1 and last_profit <= 0 and second_last_profit <= 0) {
                periods--;
                last_profit = second_last_profit;
                second_last_profit = book->revenue(created+periods-2) - cost_market - cost_unit*book->sales(created+periods-2);
            }

            // If there weren't at least 2 meaningful sales periods, we can't use this book to infer
            // anything.
            if (periods <= 1) continue;

            // TODO: this could potentially change to create new profit stream age models as we
            // encounter books of longer ages.  One tricky thing to consider: if we have a book of
            // age 6, should we really use, say, a n=1, age=6 model when we have, say, n=10, age=5
            // and n=100, age=4 models also available?  Choose one?  Weighted average?
            for (auto a : profit_stream_ages) {
                if (periods > a)
                    to_learn[a].push_back(book);
                else
                    break;
            }
        }
    }

    double weaken = creativity_.priorWeight();
    if (weaken != 1.0) {
        // Weaken existing PS beliefs
        for (auto &psb : profit_stream_beliefs_) {
            psb.second = std::move(psb.second).weaken(weaken);
        }
    }

    for (auto &learn : to_learn) {
        const auto &age = learn.first;
        RowVectorXd y(learn.second.size());
        MatrixXd X(learn.second.size(), age);

        size_t row = 0;
        for (auto &book : learn.second) {
            double total_profit = 0;
            // Figure out total profit for the book, but skip any trailing negative profit periods
            // In other words, figure out profits under the optimal removal decision.
            eris_time_t last_profit_t = 0; // Will be > 0 once we've figured it out
            for (auto t = book->outOfPrint()-1; t >= book->created(); t--) {
                double prof_t = book->revenue(t) - cost_market - cost_unit * book->sales(t);
                if (last_profit_t == 0) {
                    if (prof_t <= 0) continue;
                    // Else this book had positive profits in t, so stop skipping
                    last_profit_t = t;
                }
                total_profit += prof_t;
            }

            double cumul_profit = 0;
            const unsigned long &created = book->created();
            for (unsigned int i = 0; i < age; i++) {
                eris_time_t t = created + i;
                double prof_t = book->revenue(t) - cost_market - cost_unit * book->sales(t);
                X(row, i) = prof_t;
                cumul_profit += prof_t;
            }
            // Handle the special case described above where the last period is negative but the second
            // last is positive; in such a case, total_profit == cumul_profit, but we actually want
            // to use the final negative profit value instead of 0.  Otherwise, we want the
            // difference between total_profit and cumul_profit (== profit remaining).
            if (last_profit_t == created + age - 1)
                y[row] = book->revenue(created+age-1) - cost_market - cost_unit*book->sales(created+age-1);
            else
                y[row] = total_profit - cumul_profit;
            row++;
        }

        // Create a new non-informative belief of the given size if none exists yet
        if (profit_stream_beliefs_.count(age) == 0) {
            profit_stream_beliefs_.emplace(age, ProfitStream(age));
        }

        // Update the belief (existing or just-created) with the new data
        profit_stream_beliefs_.at(age) = profit_stream_beliefs_.at(age).update(y, X);
    }
*/
}

void Reader::updateProfitBelief() {
    // FIXME: need to handle private/public markets

    const eris_time_t t = simulation()->t();

    // If we aren't starting simulation period t=3 or later (plus the creation lag), don't do
    // anything: we can't incorporated books into the belief because we need the lagged market book
    // count.
    // +3 because nothing happens in 0, then we need a period with some creation, and the t has
    // already been incremented for the upcoming period.
    if (t < creativity_.parameters.creation_time + 3) return;

    std::list<SharedMember<Book>> new_prof_books, extrap_books;

    // Also need to go through books that we don't own, using predicted quality
    for (auto &b : book_cache_market_) {
        if (b->hasPrivateMarket()) {
            // The book is still on the market, so we'll have to extrapolate using profit stream
            // beliefs
            extrap_books.push_back(b);
        }
        else if (creativity_.publicSharingActive()) {
            // If public sharing exists, we wait until the book has been off the private market for
            // 5 periods (so that it has most likely earned any public money it's going to get)
            if (t - b->leftPrivateMarket() >= 5)
                // The book just left the market
                new_prof_books.push_back(b);
        }
        else {
            // No public market: learn as soon as it leaves the market
            new_prof_books.push_back(b);
        }
    }

    double weaken = creativity_.priorWeight();

    if (not new_prof_books.empty()) {
        MatrixXd X(new_prof_books.size(), profit_belief_->K());
        VectorXd y(new_prof_books.size());

        size_t i = 0;
        for (auto &book : new_prof_books) {
            // Calculate the book's total (private) profit and any public prize money
            y[i] = book->lifeProfitPrivate() + book->lifePrize();
            X.row(i) = Profit::profitRow(book->qualityMean(), book->order(), creativity_.market_books_avg);
            i++;
        }

        *profit_belief_ = Profit(std::move(*profit_belief_), y, X, weaken);
    }
    else if (weaken != 1.0) {
        *profit_belief_ = Profit(std::move(*profit_belief_), weaken);
    }

    // Also include books still on the market in an extrapolated belief, but we want to throw these
    // away later (because their profit level could change)
    if (extrap_books.empty()) {
        // No extrapolation books, so the "extrapolation" belief is just the profit belief
        profit_belief_extrap_ = profit_belief_;
    }
    else {
        MatrixXd X(extrap_books.size(), profit_belief_->K());
        VectorXd y(extrap_books.size());

        size_t i = 0;
        for (auto &book : extrap_books) {
            // Calculate the book's total profit (ignoring initial creation cost)
            y[i] = book->lifeProfitPrivate() + book->lifePrize();
            X.row(i) = Profit::profitRow(book->qualityMean(), book->order(), creativity_.market_books_avg);
            i++;
        }

        // extrapolation uses just-updated non-extrapolation as prior, with no weakening
        profit_belief_extrap_.reset(new Profit(*profit_belief_, y, X));
    }

}

void Reader::alterUtility(double amount) {
    auto stage = simulation()->runStage();
    if (stage == Simulation::RunStage::intra_Finish or (stage == Simulation::RunStage::intra_Apply and simulation()->runStagePriority() > 0)) {
        u_curr_ -= amount;
        u_lifetime_ -= amount;
    }
    else {
        throw std::logic_error("Reader::alterUtility can only be called in intra-finish or late-stage intra-apply simulation stages");
    }
}

void Reader::intraInitialize() {
    auto nb = creativity_.newBooks();
    for (auto &bm : nb.first) {
        book_cache_.insert(bm);
        if (bm->hasPrivateMarket()) book_cache_market_.insert(bm);
    }
}

void Reader::intraOptimize() {
    auto lock = readLock();

    std::vector<SharedMember<Book>> cache_del;

    // We get the list of all books that are available either from the market, or via piracy, and
    // calculate the utility we will receive from reading each book (but *not* including the utility
    // reduction that occurs from reading multiple books).  For purchasable books, we store the net
    // utility (i.e. reading utility minus market price).  For pirated books, we store the net
    // utility (utility minus piracy cost).  Note, however, that when the piracy catching policy
    // comes into effect, we actually have to adjust the latter by the expected increase in expected
    // piracy catch cost and fines, which makes the calculation trickier.

    // These are sorted by net utility in highest-to-lowest order for the purposes of iterating from
    // best to worse options.  We also use a multimap rather than a map here because it is
    // technically possible, if extremely unlikely, that two books have the same net utility value.
    std::multimap<double, SharedMember<Book>, std::greater<double>> market_books_u,
        piracy_books_u;

    double piracy_cost = piracyCost();
    for (auto &book : book_cache_) {
        if (library_.count(book) > 0) {
            // Already have the book: remove from cache and skip this book
            cache_del.push_back(book);
            continue;
        }
        double u = uBook(book);
        // Two ways to obtain the book:
        // - If it's on the market, can buy at its current price
        // - If piracy has been invented and one of my friends has it, I can get a copy from the friend

        if (book->hasAnyMarket()) { // On the market
            double p = book->market()->price();
            if (u > p)
                market_books_u.emplace(u - p, book);
        }

        if (creativity_.piracy() and u > piracy_cost) {
            // Look through all my friends to see if any have a copy of the book that I can copy
            for (auto &f : friends()) {
                // I can obtain a pirated book if my friend has it and my friend isn't the author
                // (or is the author, but has removed the book from the market).
                if (f->library().count(book) > 0 and (f != book->author() or not book->hasPrivateMarket())) {
                    piracy_books_u.emplace(u - piracy_cost, book);
                    break;
                }
            }
        }
    }

    lock.write(); // Lock out other optimizers from accessing this agent until we're done modifying things.

    for (auto &del : cache_del) book_cache_.erase(del);

    // Get expected piracy penalty value:
    SharedMember<CopyrightPolice> police;
    auto exp_piracy_cost = [this, &police](unsigned pirated) -> double {
        if (!creativity_.catchPiratesActive())
            return 0.0;

        double prob = police->prob(pirated);
        double cost = police->cost() + police->fine(pirated);
        return prob*cost;
    };

    unsigned books_pirating = 0;
    double piracy_penalty = 0, piracy_penalty_next = 0;
    if (creativity_.catchPiratesActive())
        police = simulation()->agents<CopyrightPolice>()[0];

    // Baseline expected penalty: even if we don't pirate, there's some probability of being
    // (falsely) caught for piracy.  (NB, if we don't try catching pirates, this is simply 0).
    piracy_penalty = exp_piracy_cost(0);

    // The expected cost of piracy if pirating an additional book (also 0 if not catching pirates).
    piracy_penalty_next = exp_piracy_cost(1);

    // Now we're going to iterate through the market and piracy books, taking one from whichever is
    // better in each pass of this loop.  For piracy books, we have to consider not only the net
    // utility value, but also, when the catch policy is active, the increase in the expected cost
    // and fine of being caught pirating.
    auto market_it = market_books_u.begin();
    auto piracy_it = piracy_books_u.begin();
    std::set<SharedMember<Book>> new_books;
    double money = assets[creativity_.money];
    double u_curr = u(money, new_books) - piracy_penalty; // the "no-books" utility

    while (market_it != market_books_u.end() || piracy_it != piracy_books_u.end()) {
        double market_u = -std::numeric_limits<double>::infinity(), piracy_u = -std::numeric_limits<double>::infinity();

        // If we've already acquired this book (i.e. because the book is available both via piracy
        // and the market, and we already obtained one and are now encountering the other), skip.
        while (market_it != market_books_u.end() and new_books.count(market_it->second))
            ++market_it;
        while (piracy_it != piracy_books_u.end() and new_books.count(piracy_it->second))
            ++piracy_it;

        if (market_it != market_books_u.end())
            market_u = market_it->first;
        if (piracy_it != piracy_books_u.end())
            piracy_u = piracy_it->first - (piracy_penalty_next - piracy_penalty);

        SharedMember<Book> book;
        double cost; // Direct cost: either the market price, or the piracy cost
        bool piracy = false;
        // Now consider the best market book with the best piracy book: we'll consider whichever
        // gives us the higher utility increase (including the potential increase in expected fine,
        // for the piracy option).
        if (market_u > 0 and market_u >= piracy_u) { // The best market option is better than the best piracy option
            book = market_it->second;
            cost = book->market()->price();
            ++market_it;
        }
        else if (piracy_u > 0) { // Otherwise the best piracy option gives the highest utility
            book = piracy_it->second;
            cost = piracy_cost;
            piracy = true;
            ++piracy_it;
        }
        else {
            // The best choice doesn't increase utility, so no point in looking further (later books
            // have a net utility no higher than this one).  Note that we don't typically get here:
            // this 'else' is only hit if piracy and pirate catching are enabled, and there are
            // positive utility piracy books, but the increase in penalty exceeds the utility gain
            // of the book.  (We never add on-market books with negative utility in the first place,
            // and piracy books are only added when the utility exceeds the direct piracy cost).
            break;
        }

        if (cost > money) {
            // We can't actually afford the direct cost of the book (i.e. the market price, or the
            // direct piracy cost), so skip it but don't stop looking: there may be less good books
            // that we can still afford.  For example, imagine a book that gives 1003 utils but
            // costs 1000, which will be considered before a book that gives 3 utils and costs 1.
            //
            // Note also that running out of income is potentially a misoptimization: in the above
            // example, if the individual has income of 1000, and there are 1000 3 util/1 cost books
            // available, with this optimization he'll buy the single 1000-cost book and get a net
            // utility gain of 3, but if he bought the 1000 1-cost books he'd had net utility gain
            // of 2000.  In practice, this limitation is not particularly relevant as income is
            // typically made high enough to not often bind.
            continue;
        }

        new_books.insert(book);

        // Get expected utility with the new book included (i.e. known utility minus the expected
        // piracy penalty).  This isn't quite the same as the above because there may be an
        // additional opportunity cost for adding another book (i.e. a too-many-books penalty).
        double u_with_book = u(money - cost, new_books) - (piracy ? piracy_penalty_next : piracy_penalty);

        if (u_with_book <= u_curr) {
            // Obtaining this book doesn't increase utility, so skip it.
            new_books.erase(book);
            continue;
        }

        // Otherwise the utility with the book, after giving up the cost of the book, the
        // opportunity cost of reading an additional book, and the expected cost of getting caught
        // pirating with an extra book, is higher than without the book, so go ahead with the
        // purchase/pirating.

        u_curr = u_with_book;
        money -= cost;

        if (!piracy) { // Buying from the market
            auto bm = book->market();
            auto bm_lock = lock.supplement(bm);

            reservations_.push_front(bm->reserve(sharedSelf(), 1, cost));
            reserved_books_.emplace(book, bm->isPublic() ? BookCopy::Status::purchased_public : BookCopy::Status::purchased_market);
        }
        else { // Pirating
            books_pirating++;
            piracy_penalty = piracy_penalty_next;
            piracy_penalty_next = exp_piracy_cost(books_pirating + 1);
            reserved_books_.emplace(book, BookCopy::Status::pirated);
            reserved_piracy_cost_ += piracy_cost;
        }
    }
}

void Reader::intraApply() {
    auto lock = writeLock();
    // Apply the reserved transactions (i.e. pay for market book copies and receive the book)
    for (auto &res : reservations_) {
        res.buy();
    }
    // Also remove any cost associated with obtaining pirated copies of books:
    assets.transferApprox({ creativity_.money, reserved_piracy_cost_ }, 1e-6);
    reserved_piracy_cost_ = 0.0;
    reservations_.clear();

    // Reset the "new" books list to the set of just-bought and just-pirated books
    library_new_.clear();

    for (auto &new_book : reserved_books_) {
        auto &book = new_book.first;
        auto &status = new_book.second;

        // Remove the book from assets (where the reservation puts it) -- we track it via library_ instead.
        assets.remove(book);

        auto inserted = library_.emplace(
                SharedMember<Book>(book),
                BookCopy(book->qualityDraw(), status, simulation()->t()));

        // And it's always considered a "new" (to us) book, so stash it there:
        library_new_.emplace(book, std::ref(inserted.first->second));

        // If we obtained a pirated book, record that in the base book metadata:
        if (status == BookCopy::Status::pirated) {
            auto lock_b = lock.supplement(book);
            book->recordPiracy(1);
        }
        //else {} -- // book->recordSale() not needed here: it's called during the BookMarket transaction
    }
    reserved_books_.clear();

    if (creativity_.publicVotingActive()) {
        // If we are voting, figure out the distribution of votes, by allocating votes in coarse
        // proportion to the realized book sub-utility (i.e. uBook(), the book utility, minus its
        // cost (market or piracy), without taking into account the added opportunity cost of
        // reading additional books).
        //
        // To do that, we first order all public-sharing-purchased books by their sub-utility value
        // (highest first).  Then we assign a vote to the first and reweight its subutility.  We
        // repeat this until all votes are used up.
        //
        // The reweighting follows the pattern 1/2, 2/3, 3/4, ... (compounded, so that the stored
        // subutility value is 1/2, 1/3, 1/4, ... of the original subutility).  Effectively this
        // works by ensuring that we don't give a second vote to a book unless its utility is at
        // least double the next-best-book's utility, and hold back a third vote unless it is at
        // least three times the next-best-book's utility.
        //
        // For example, with book/utility pairs: (A,10), (B,4), (C,3), we would assign votes in this
        // order (stopping when we run out of votes; this pattern repeats indefinitely):
        // A, A, B, A, C, A, [AB], A, C, A, B, A, A, [ABC]
        // (where [XYZ] indicates a sequence with any subordering of X, Y, Z).
        //
        // Note that, asymptotically, this tends to award a vote distribution equal to the relative
        // utility distribution (in this special case of integer utilities, matching exactly after
        // 17n votes).

        uint32_t total_votes = 0;
        std::unordered_map<SharedMember<Book>, uint32_t> votes;
        std::multimap<double, SharedMember<Book>, std::greater<double>> book_subu;

        for (auto &bc : library_new_) {
            auto &copy = bc.second.get();
            if (copy.purchased_public()) {
                auto book = bc.first;
                double netu = uBook(book, copy.quality) - book->price();
                if (std::isfinite(netu) and netu > 0)
                    book_subu.emplace(std::move(netu), std::move(book));
            }
        }

        if (not book_subu.empty()) {

            while (total_votes < creativity_.parameters.policy_public_voting_votes) {
                // We need to change the key, which means we have to remove and reinsert the element
                // (so that it gets resorted into the right place in the multimap).
                auto best = book_subu.begin();
                auto book = best->second;
                auto &v = votes[book];
                v++;
                total_votes++;
                double netu_updated = (best->first * v) / (v+1);
                book_subu.erase(best);
                book_subu.emplace(std::move(netu_updated), std::move(book));
            }

            for (auto &bv : votes) {
                auto lock_b = lock.supplement(bv.first);
                bv.first->vote(bv.second);
            }
        }
    }

    // "Eat" any money leftover
    double money = assets.remove(creativity_.money);
    // Store final utility
    u_curr_ = u(money, library_new_);
    u_lifetime_ += u_curr_;
}

void Reader::intraReset() {
    auto lock = writeLock();
    reservations_.clear();
    reserved_books_.clear();
    reserved_piracy_cost_ = 0;
}

}
