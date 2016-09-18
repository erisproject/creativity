#include "creativity/Creativity.hpp"
#include "creativity/Reader.hpp"
#include "creativity/Book.hpp"
#include "creativity/BookMarket.hpp"
#include "creativity/belief/Profit.hpp"
#include <eris/debug.hpp>
#include <eris/random/util.hpp>
#include <boost/random/chi_squared_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
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
    double plast = evalPolynomial(0, coef);
    for (const double &x : {1e-100, 1e-10, .1, 1., 10., 100., 1e10}) {
        double p = evalPolynomial(x, coef);
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
    return evalPolynomial(distance, distancePenaltyPolynomial());
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
    double plast = evalPolynomial(xlast, coef);
    for (const unsigned &x : {1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,1000,1000000}) {
        double p = evalPolynomial(x, coef);
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
    return evalPolynomial(n, numBooksPenaltyPolynomial());
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
#   ifdef ERIS_DEBUG
    try {
#   endif
    // Update the various beliefs
    updateBeliefs();

    double income_available = assets[creativity_.money] + creativity_.parameters.income;
    if (creativity_.publicSharing()) income_available -= creativity_.parameters.policy_public_sharing_tax;
    if (creativity_.catchPirates()) income_available -= creativity_.parameters.policy_catch_tax;

    auto sim = simulation();
    const double market_books = sim->countMarkets<BookMarket>();
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
                            sim->t() - 1 - book->lastSale(), book->age(), authored_books - 1, market_books, cost_unit,
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
        if (creativity_.parameters.creation_time == 0) creation_base_cost += cost_market;

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
        if (creativity_.publicSharing() and creativity_.parameters.cost_piracy < marginal_cost)
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
                            authored_books, market_books,
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
                            creativity_.parameters.readers, 0, 0, authored_books, market_books,
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
        }
    }

    if (create_countdown_ == 0) {
        if (creativity_.parameters.creation_time > 0) {
            // A currently-in-progress book is done (and wasn't created instantaneously), so choose its
            // price based on current demand beliefs.  If the optimal price is 0, we won't release it at
            // all: the creation is a sunk cost: it apparently doesn't seem profitable any more.  If we
            // can't draw a price at all, we'll just fall back on the initial one decided when the book
            // was created.
            if (usableBelief(demand_belief_)) {
                // As a safety against absurd prices, we also enforce a price ceiling of 20% of
                // per-period income.
                std::pair<double, double> max;
                try {
                    max = demand_belief_.argmaxP(creativity_.parameters.prediction_draws, create_quality_, 0, 0,
                            creativity_.parameters.readers, 0, 0, authored_books, market_books,
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
}

void Reader::interApply() {
    // Give potential income (this has to be before authorship decision, since authorship requires
    // giving up some potential income: actual income will be this amount minus whatever is given up
    // to author)
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
            X.row(i) = Demand::bookRow(b, b->qualityMean(), creativity_.market_books_lagged);
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
        else if (creativity_.publicSharing()) {
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
            X.row(i) = Profit::profitRow(book->qualityMean(), book->order(), creativity_.market_books_lagged);
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
            X.row(i) = Profit::profitRow(book->qualityMean(), book->order(), creativity_.market_books_lagged);
            i++;
        }

        // extrapolation uses just-updated non-extrapolation as prior, with no weakening
        profit_belief_extrap_.reset(new Profit(*profit_belief_, y, X));
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
    // Map utility-minus-price values to `<buy, book>` pairs, where `buy` is true if the book is to
    // be purchased and false if the book is to be pirated from a friend (we don't store *which*
    // friend, but it'll only be true if there is at least one friend with a copy).
    //
    // This is sorted by net utility in highest-to-lowest order for the purposes of iterating from
    // best to worse options.
    std::map<double, std::vector<std::pair<bool, SharedMember<Book>>>, std::greater<double>> book_net_u;
    double piracy_cost = piracyCost();
    for (auto &book : book_cache_) {
        if (library_.count(book) > 0) {
            // Already have the book: remove from cache and continue
            cache_del.push_back(book);
            continue;
        }
        double u = uBook(book);
        // Two ways to obtain the book:
        // - If it's on the market, can buy at its current price
        // - If piracy has been invented and one of my friends has it, I can get a copy from the friend

        bool for_sale = book->hasAnyMarket();
        bool shared = false;
        for (auto &f : friends()) {
            if (f->library().count(book) > 0) {
                // My friend owns the book
                if (book->author() != f or not book->hasPrivateMarket()) {
                    // ... and either isn't the author, or is (and the book has left the market)
                    shared = true;
                    break;
                }
            }
        }

        double p = for_sale ? book->market()->price() : 0.0;

        if (for_sale and shared and p > piracy_cost)
            for_sale = false; // It's for sale, but the pirated version is cheaper.

        if (for_sale and u >= p) // We're buying, and the book is worth buying
            book_net_u[u - p].emplace_back(true, book);
        else if (shared and u >= piracy_cost) // We're pirating, and the book is worth the cost of pirating
            book_net_u[u - piracy_cost].emplace_back(false, book);
    }

    // Shuffle the order of any books with the same net utilities so that we'll be choosing randomly
    // from equally-preferred options.
    for (auto &bnu : book_net_u) {
        if (bnu.second.size() > 1)
            std::shuffle(bnu.second.begin(), bnu.second.end(), random::rng());
    }

    lock.write(); // Time to get serious.

    for (auto &del : cache_del) book_cache_.erase(del);

    std::set<SharedMember<Book>> new_books;
    double money = assets[creativity_.money];
    double u_curr = u(money, new_books); // the "no-books" utility
    // Keep looking at books as long as the net utility from buying the next best book exceeds the
    // penalty that will be incurred.
#   ifdef ERIS_DEBUG
    try {
#   endif
    while (not book_net_u.empty() and money > 0) {
        // Pull off the best remaining book:
        auto &s = book_net_u.begin()->second;
        auto buy_book = std::move(s.back());
        bool buying = buy_book.first;
        auto &book = buy_book.second;
        s.pop_back();
        // If we just pulled off the last book at the current net utility level, delete the level.
        if (s.empty()) book_net_u.erase(book_net_u.begin());

        if (buying) {
            auto bm = book->market();
            lock.add(bm);
            // Perform some safety checks:
            // - if we've already added the book to the set of books to buy, don't add it again. (This
            //   shouldn't be possible)
            if (new_books.count(book) > 0) throw std::logic_error("Internal Reader error: attempt to buy book that is already obtained!");
            // - if the BookMarket says the book is no longer feasible, ignore it. (This shouldn't
            //   happen, but just in case).
            auto pinfo = bm->price(1);
            if (not pinfo.feasible) throw std::logic_error("Internal Reader error: book market says book not available!");
            if (pinfo.total <= money) { // Only buy if we can actually afford it
                new_books.insert(book);
                double u_with_book = u(money - pinfo.total, new_books);
                if (u_with_book < u_curr) {
                    // The best book *lowered* utility (because of the multiple books penalty), so
                    // we're done (because all other new books have utility no higher than this
                    // one), and will incur the same utility penalty.
                    new_books.erase(book);
                    lock.remove(bm);
                    break;
                }
                // Otherwise the purchase is a good one, so make a reservation.
                reservations_.push_front(bm->reserve(sharedSelf(), 1, pinfo.total));
                reserved_books_.emplace(book, bm->isPublic() ? BookCopy::Status::purchased_public : BookCopy::Status::purchased_market);
                u_curr = u_with_book;
                money -= pinfo.total;
            }
            lock.remove(bm);
        }
        else {
            // Pirating.
            if (new_books.count(book) != 0) throw std::logic_error("Internal Reader error: attempt to pirate book that is already obtained!");

            if (piracy_cost <= money) { // Only pirate if we can actually afford piracy
                new_books.insert(book);

                double u_with_book = u(money - piracy_cost, new_books);
                if (u_with_book < u_curr) {
                    // The best book *lowered* utility (because of the multiple books penalty), so
                    // we're done (because all other new books have utility no higher than this
                    // one), and will incur the same utility penalty.
                    new_books.erase(book);
                    break;
                }
                reserved_books_.emplace(book, BookCopy::Status::pirated);
                reserved_piracy_cost_ += piracy_cost;
                u_curr = u_with_book;
                money -= piracy_cost;
            }
        }
    }
#   ifdef ERIS_DEBUG
    } catch (std::exception &e) {
        std::cerr << "Caught exception in reader intraOptimize buy loop: " << e.what();
        throw e;
    }
#   endif
}
void Reader::intraApply() {
    auto lock = writeLock();
    // Apply the reserved transactions (i.e. pay for book copies and receive the book)
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

        // Remove the book from assets -- we track it via library_ instead.
        assets.remove(book);

        auto inserted = library_.emplace(
                SharedMember<Book>(book),
                BookCopy(book->qualityDraw(), status, simulation()->t()));

        // And it's always considered a "new" (to us) book, so stash it there:
        library_new_.emplace(book, std::ref(inserted.first->second));

        // If we obtained a pirated book, record that in the base book metadata:
        if (status == BookCopy::Status::pirated) book->recordPiracy(1);
        //else {} -- // book->recordSale() not needed here: it's called during the BookMarket transaction
    }
    reserved_books_.clear();

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
