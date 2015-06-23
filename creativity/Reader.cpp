#include "creativity/Creativity.hpp"
#include "creativity/Reader.hpp"
#include "creativity/Book.hpp"
#include "creativity/BookMarket.hpp"
#include <eris/Random.hpp>
#include <algorithm>
#include <map>
#include <list>
#include <random>
#include <cmath>

using namespace eris;
using namespace Eigen;
using namespace eris::belief;
using namespace creativity::belief;

namespace creativity {

const std::vector<double> Reader::default_polynomial{{4., -1.}};
const std::vector<double> Reader::default_penalty_polynomial{{0, 0, 0.25}};
const std::vector<unsigned int> Reader::profit_stream_ages{{1,2,4,8}};

Reader::Reader(
        std::shared_ptr<Creativity> creativity,
        const Position &pos,
        double cFixed, double cUnit, double cPiracy, double inc
        )
    : WrappedPositional<agent::AssetAgent>(pos, creativity->parameters.boundary, -creativity->parameters.boundary),
    cost_fixed{cFixed}, cost_unit{cUnit}, cost_piracy{cPiracy}, income{inc},
    creativity_{std::move(creativity)},
    profit_belief_{new Profit()}, profit_belief_extrap_{profit_belief_}
{
    profit_stream_beliefs_.emplace(1, ProfitStream(1));
}

void Reader::uPolynomial(std::vector<double> coef) {
    // 0th coefficient must be positive
    if (coef.size() >= 1 and coef[0] <= 0) throw std::domain_error("Invalid uPolynomial: coef[0] <= 0");

    // all coefficients must be finite; first and last non-zero coefficients must be negative
    double first = 0.0, last = 0.0;
    for (size_t i = 0; i < coef.size(); i++) {
        if (not std::isfinite(coef[i]))
            throw std::domain_error("Invalid uPolynomial: coef[" + std::to_string(i) + "] is not finite");

        if (i > 0 and coef[i] != 0) {
            last = coef[i];
            if (first == 0) first = coef[i];
        }
    }

    if (first > 0) throw std::domain_error("Invalid uPolynomial: first non-zero, non-constant coefficient must be negative.");
    if (last > 0) throw std::domain_error("Invalid uPolynomial: last non-zero, non-constant coefficient must be negative.");

    double xlast = 0;
    double p = evalPolynomial(0, coef);
    for (const double &x : {1e-6, .001, .01, .1, 1., 10., 100., 1000., 1e6}) {
        double pnext = evalPolynomial(x, coef);
        if (pnext > p)
            throw std::domain_error("Invalid uPolynomial: polynomial is increasing: f(" + std::to_string(x) + ") > f(" + std::to_string(xlast) + ")");
    }

    u_poly_ = std::move(coef);
}

double Reader::u(double money, const std::unordered_map<eris::SharedMember<Book>, std::reference_wrapper<BookCopy>> &books) const {
    double u = money - penalty(books.size());
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

const std::vector<double>& Reader::uPolynomial() const {
    return u_poly_;
}

double Reader::evalPolynomial(double x, const std::vector<double> &polynomial) {
    double p = 0.0;
    double xi = 1.0;
    for (auto &c : polynomial) {
        p += c * xi;
        xi *= x;
    }
    return p;
}

void Reader::penaltyPolynomial(std::vector<double> coef) {
    double last = evalPolynomial(0, coef);
    unsigned long last_x = 0;
    for (auto &x : {1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,1000,1000000}) {
        double pen = evalPolynomial(x, coef);
        if (pen < last)
            throw std::domain_error("Invalid penalty polynomial: f(" + std::to_string(x) + ") < f("
                    + std::to_string(last_x) + ")");
        last = pen;
        last_x = x;
    }

    pen_poly_ = std::move(coef);
}
const std::vector<double>& Reader::penaltyPolynomial() const {
    return pen_poly_;
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
    return creativity_->parameters.cost_piracy;
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

void Reader::interOptimize() {
#   ifdef ERIS_DEBUG
    try {
#   endif
    // Update the various profit, demand, and quality beliefs
    updateBeliefs();

    double income_available = assets()[creativity_->money] + income;

    auto sim = simulation();
    const double market_books = sim->countMarkets<BookMarket>();
    const double authored_books = wrote().size();

    // Any books not in new_prices_ but that were on the market last period will be removed from the
    // market
    new_prices_.clear();

    auto &rng = Random::rng();

    std::list<eris::SharedMember<Book>> bypass_beliefs;
    if (usableBelief(demand_belief_)) {
        // Figure out the expected profitability of each book; positive ones will be kept on the market
        std::map<double, std::vector<std::pair<SharedMember<Book>, double>>, std::greater<double>> profitability;
        for (auto &book : wrote_market_) {
            auto book_mkt = book->market();

            try {
                auto max = demand_belief_.argmaxP(creativity_->parameters.prediction_draws, book->quality(), book->lifeSales(), book->lifePirated(), creativity_->parameters.readers,
                        sim->t() - 1 - book->lastSale(), book->age(), authored_books - 1, market_books, cost_unit,
                        income / 10);
                const double &p = max.first;
                const double &q = max.second;

                const double profit = (p - cost_unit) * q - cost_fixed;

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
        while (income_available >= cost_fixed and not profitability.empty()) {
            auto &books = profitability.begin()->second;
            if (income_available < cost_fixed * books.size())
                // We don't have enough to do all the books at this profitability level, so shuffle them
                // so that we choose randomly
                std::shuffle(books.begin(), books.end(), rng);

            for (auto &b : books) {
                if (income_available < cost_fixed) break;
                income_available -= cost_fixed;
                new_prices_.emplace(b.first, b.second);
            }
        }
    }
    else {
        bypass_beliefs.insert(bypass_beliefs.end(), wrote_market_.begin(), wrote_market_.end());
    }

    // Bypass market predictions if our market demand belief is noninformative, not sufficiently
    // informed, or can't produce usable draws.  In these cases, use initial model action
    // parameters instead.
    if (creativity_->parameters.initial.prob_keep > 0) {
        std::bernoulli_distribution keep(creativity_->parameters.initial.prob_keep);
        for (auto on_market = bypass_beliefs.cbegin();
                income_available >= cost_fixed and on_market != bypass_beliefs.cend(); on_market++) {
            auto &book = *on_market;
            if (keep(rng)) {
                income_available -= cost_fixed;
                double new_price = (book->price() - cost_unit) * creativity_->parameters.initial.keep_price + cost_unit;
                new_prices_.emplace(book, new_price);
            }
        }
    }

    if (create_countdown_ == -1) {
        //// NEW BOOK CREATION:
        if (income_available >= cost_fixed) {
            // Creating a book requires an ability to predict profit (to determine whether creation is
            // worthwhile), and an ability to predict demand (to determine the initial price)
            if (usableBelief(*profit_belief_) and usableBelief(demand_belief_)) { // NB: profit_belief_extrap_ is a copy or update of profit_belief_, so will also be usable
                // Create a book if E(profit) > effort required

                // Find the l that maximizes creation profit given current wealth (which would have come from
                // past sales, if non-zero) plus the fixed income we're about to receive, minus whatever we
                // decided to spend above to keep books on the market.

                try {
                    double effort;
#               ifdef ERIS_DEBUG
                    try {
#               endif
                    effort = profit_belief_extrap_->argmaxL(
                            creativity_->parameters.prediction_draws,
                            [this] (double l) -> double { return creationQuality(l); },
                            authored_books, market_books, income_available - cost_fixed
                            );
#               ifdef ERIS_DEBUG
                    } catch (BayesianLinear::draw_failure &e) {
                        ERIS_DBG("draw failure in profit_extrap argmaxL for reader=" << id() << ", t=" << simulation()->t() << ": " << e.what());
                        throw;
                    }
#               endif
                    double quality = creationQuality(effort);

                    double exp_profit;
#               ifdef ERIS_DEBUG
                    try {
#               endif
                    exp_profit = profit_belief_extrap_->predict(creativity_->parameters.prediction_draws, quality, authored_books, market_books);
#               ifdef ERIS_DEBUG
                    }
                    catch (BayesianLinear::draw_failure &e) {
                        ERIS_DBG("draw failure in profit_extrap belief prediction; reader="<<id() << ", t=" << simulation()->t());
                        throw;
                    }
#               endif

                    // If the optimal effort level gives a book with positive expected profits, do it:
                    if (exp_profit > effort) {
                        // We're going to create, so calculate the optimal first-period price.  It's possible that
                        // we get back cost_unit (and so predicted profit is non-positive); write the book anyway:
                        // perhaps profits are expected to come in later periods?
                        std::pair<double, double> max;
#                   ifdef ERIS_DEBUG
                        try {
#                   endif
                        max = demand_belief_.argmaxP(creativity_->parameters.prediction_draws, quality, 0, 0,
                                creativity_->parameters.readers, 0, 0, authored_books, market_books, cost_unit,
                                income / 10);
#                   ifdef ERIS_DEBUG
                        }
                        catch (BayesianLinear::draw_failure &e) {
                            ERIS_DBG("draw failure in demand argmaxP calculation; reader="<<id() << ", t=" << simulation()->t());
                            ERIS_DBGVAR(demand_belief_.draw_rejection_success);
                            throw;
                        }
#                   endif
                        if (max.first > 0) {
                            create_price_ = max.first;
                            create_quality_ = quality;
                            create_effort_ = effort;
                            create_started_ = true;
                            create_countdown_ = creativity_->parameters.creation_time;
                            create_position_ = position();
                        }
                        // else: even though our Profit seems like it would be positive, our demand
                        // belief suggests otherwise, so don't create.
                    }
                }
                catch (BayesianLinear::draw_failure &e) {
                    // Ignore draw failures
                }
            }
            else {
                // If we have no useful profit belief yet, just use the initial values:
                if (creativity_->parameters.initial.prob_write > 0 and std::bernoulli_distribution(creativity_->parameters.initial.prob_write)(rng)) {
                    double q = std::uniform_real_distribution<double>(creativity_->parameters.initial.q_min, creativity_->parameters.initial.q_max)(rng);
                    double effort = creationEffort(q);
                    // Make sure the required effort doesn't exceed the available funds
                    if (income_available >= effort) {
                        create_started_ = true;
                        create_countdown_ = creativity_->parameters.creation_time;
                        create_price_ = cost_unit + std::uniform_real_distribution<double>(creativity_->parameters.initial.p_min, creativity_->parameters.initial.p_max)(rng);
                        create_quality_ = q;
                        create_effort_ = effort;
                        create_position_ = position();
                    }
                }
            }
        }
    }
#   ifdef ERIS_DEBUG
    } catch (std::exception &e) {
        std::cerr << "\nCaught exception in reader interOptimize: " << e.what() << "\n";
        throw;
    }
#   endif
}

void Reader::interApply() {
    // Give potential income (this has to be before authorship decision, since authorship requires
    // giving up some potential income: actual income will be this amount minus whatever is given up
    // to author)
    assets()[creativity_->money] += income;

    auto sim = simulation();
    SharedMember<Book> newbook;
    if (create_started_) {
        // Incur the creation effort cost right away
        assets()[creativity_->money] -= create_effort_;
        create_started_ = false;
    }
    if (create_countdown_ > 0) {
        create_countdown_--;
    }
    else if (create_countdown_ == 0) {
        // Book is finished, release it.

        if (assets()[creativity_->money] >= cost_fixed) {
            // We can't afford to release it right now, so just hang onto it and release it next
            // period.

            // The cost (think of this as an opportunity cost) of creating, and the first period fixed
            // cost:
            assets()[creativity_->money] -= cost_fixed;

            auto qdraw = [this] (const Book &book) -> double {
                // Truncated normal.  Since book.quality() > 0, this should return a valid draw 50% of
                // the time, so this loop shouldn't run that much, usually.
                double x;
                do { x = book.quality() + writer_quality_sd * Random::rstdnorm(); }
                while (x < 0);
                return x;
            };
            newbook = sim->spawn<Book>(creativity_, create_position_, sharedSelf(), wrote_.size(), create_price_, create_quality_, qdraw);

            /// If enabled, add some noise in a random direction to the position
            if (writer_book_sd > 0) {
                double step_dist = writer_book_sd * Random::rstdnorm();
                newbook->moveBy(step_dist * Position::random(create_position_.dimensions));
            }

            wrote_.insert(wrote_.end(), newbook);
            wrote_market_.insert(newbook);
            auto ins = library_.emplace(SharedMember<Book>(newbook), BookCopy(create_quality_, BookCopy::Status::wrote, sim->t()));
            library_on_market_.emplace(newbook, std::ref(ins.first->second));

            create_countdown_--;
        }
        // Otherwise we can't afford to bring it to market just now, so hold onto it until next
        // time.
    }

    std::vector<SharedMember<Book>> remove;
    for (auto &b : wrote_market_) {
        // If it's staying on the market, update the price and incur the fixed cost
        if (new_prices_.count(b) > 0) {
            b->market()->setPrice(new_prices_[b]);
            assets()[creativity_->money] -= cost_fixed;
        }
        else if (newbook != b) {
            // No new price (and not the newbie), which means we're removing the book from the market
            remove.push_back(b);
        }
    }

    for (auto &b : remove) {
        wrote_market_.erase(b);
        simulation()->remove(b->market());
    }

    // Finally, move a random distance in a random direction
    if (creativity_->parameters.reader_step_sd > 0) {
        double step_dist = Random::rstdnorm() * creativity_->parameters.reader_step_sd;
        if (step_dist != 0)
            moveBy(step_dist * Position::random(position().dimensions));
    }
}

double Reader::uBook(const SharedMember<Book> &b) const {
    return uBook(b, quality(b));
}

double Reader::uBook(const SharedMember<Book> &b, double quality) const {
    quality += evalPolynomial(distance(b), u_poly_);
    if (quality < 0) quality = 0.0;
    return quality;
}

double Reader::quality(const SharedMember<Book> &b) const {
    auto found = library_.find(b);
    if (found != library_.end())
        return found->second.quality;

    auto found_pred = quality_predictions_.find(b);
    if (found_pred != quality_predictions_.end())
        return found_pred->second;

    double q_hat;
    if (usableBelief(quality_belief_)) {
        // Use the quality belief to predict the quality
        auto &q_b = const_cast<belief::Quality&>(quality_belief_);
        q_hat = q_b.predict(b, creativity_->parameters.prediction_draws);
    }
    else {
        // No informative beliefs above quality, so just use initial parameter mean (= midpoint)
        q_hat = (creativity_->parameters.initial.q_max + creativity_->parameters.initial.q_min) / 2.0;
    }
    if (q_hat < 0) q_hat = 0;
    quality_predictions_.emplace(b, q_hat);

    return q_hat;
}

double Reader::penalty(unsigned long n) const {
    return evalPolynomial(n, pen_poly_);
}

const std::unordered_map<SharedMember<Book>, BookCopy>& Reader::library() const { return library_; }

void Reader::receiveProceeds(const SharedMember<Book> &book, const Bundle &revenue) {
    assets() += revenue;
    Bundle tvc(creativity_->money, book->currSales() * cost_unit);
    tvc.transferApprox(tvc, assets(), 1e-8);
}

const belief::Profit& Reader::profitBelief() const { return *profit_belief_; }
const belief::Profit& Reader::profitExtrapBelief() const { return *profit_belief_extrap_; }
bool Reader::profitExtrapBeliefDiffers() const { return profit_belief_ != profit_belief_extrap_; }
const belief::Demand& Reader::demandBelief() const { return demand_belief_; }
const belief::Quality& Reader::qualityBelief() const { return quality_belief_; }
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
    return not model.noninformative() and model.n() - model.K() >= creativity_->parameters.initial.belief_threshold;
}

void Reader::updateBeliefs() {
    updateQualityBelief();
    updateDemandBelief();
    //updateProfitStreamBelief();
    updateProfitBelief();

    // Clear any "on-market" book references that aren't on the market anymore, since they just got
    // incorporated into the beliefs above.
    std::vector<SharedMember<Book>> remove;
    for (auto &book : library_on_market_) {
        if (not book.first->hasMarket())
            remove.push_back(book.first);
    }
    for (auto &book : remove)
        library_on_market_.erase(book);

    remove.clear();
    for (auto &book : book_cache_market_) {
        if (not book->hasMarket())
            remove.push_back(book);
    }
    for (auto &book : remove)
        book_cache_market_.erase(book);
}

void Reader::updateQualityBelief() {
    double weaken = creativity_->priorWeight();
    bool belief_changed = true;

    // If we obtained any new books, update the quality belief with them
    if (not newBooks().empty()) {
        std::vector<SharedMember<Book>> books;
        VectorXd y(newBooks().size());
        size_t i = 0;
        for (auto &a : newBooks()) {
            books.push_back(a.first);
            y[i++] = a.second.get().quality;
        }

        quality_belief_ = std::move(quality_belief_).weaken(weaken).update(y, Quality::bookData(books));
    }
    else if (weaken != 1.0) {
        quality_belief_ = std::move(quality_belief_).weaken(weaken);
    }
    else {
        belief_changed = false;
    }

    if (belief_changed) {
        quality_predictions_.clear();
    }
}
void Reader::updateDemandBelief() {
    // If we aren't starting simulation period t=3 or later, don't do anything: we can't
    // incorporated books into the belief because we need the lagged market book count
    if (simulation()->t() < 3) return;

    double weaken = creativity_->priorWeight();

    // Figure out which books might yield usable data: i.e. those on the market
    std::list<SharedMember<Book>> mktbooks;
    for (auto &bu : library_on_market_) {
        if (bu.first->hasMarket()) mktbooks.push_back(bu.first);
    }
    for (auto &b : book_cache_market_) {
        if (b->hasMarket()) mktbooks.push_back(b);
    }

    if (not mktbooks.empty()) {
        VectorXd y(mktbooks.size());
        MatrixXdR X(mktbooks.size(), demand_belief_.K());
        // NB: this runs in the interoptimizer, which means t has already been incremented
        auto last_t = simulation()->t() - 1;
        size_t i = 0;
        for (const auto &b : mktbooks) {
            y[i] = b->sales(last_t);
            X.row(i) = Demand::bookRow(b, quality(b), creativity_->market_books_lagged);
            i++;
        }

        demand_belief_ = std::move(demand_belief_).weaken(weaken).update(y.head(i), X.topRows(i));
    }
    else if (weaken != 1.0)
        demand_belief_ = std::move(demand_belief_).weaken(weaken);
}
void Reader::updateProfitStreamBelief() {

    // Figure out which books might yield usable data
    std::list<SharedMember<Book>> potential;
    for (auto &bu : library_on_market_) {
        if (not bu.first->hasMarket()) potential.push_back(bu.first);
    }
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
            double last_profit = book->revenue(created+periods-1) - cost_fixed - cost_unit*book->sales(created+periods-1);
            double second_last_profit = book->revenue(created+periods-2) - cost_fixed - cost_unit*book->sales(created+periods-2);
            while (periods > 1 and last_profit <= 0 and second_last_profit <= 0) {
                periods--;
                last_profit = second_last_profit;
                second_last_profit = book->revenue(created+periods-2) - cost_fixed - cost_unit*book->sales(created+periods-2);
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

    double weaken = creativity_->priorWeight();
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
                double prof_t = book->revenue(t) - cost_fixed - cost_unit * book->sales(t);
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
                double prof_t = book->revenue(t) - cost_fixed - cost_unit * book->sales(t);
                X(row, i) = prof_t;
                cumul_profit += prof_t;
            }
            // Handle the special case described above where the last period is negative but the second
            // last is positive; in such a case, total_profit == cumul_profit, but we actually want
            // to use the final negative profit value instead of 0.  Otherwise, we want the
            // difference between total_profit and cumul_profit (== profit remaining).
            if (last_profit_t == created + age - 1)
                y[row] = book->revenue(created+age-1) - cost_fixed - cost_unit*book->sales(created+age-1);
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
}

void Reader::updateProfitBelief() {
    // If we aren't starting simulation period t=3 or later, don't do anything: we can't
    // incorporated books into the belief because we need the lagged market book count
    if (simulation()->t() < 3) return;

    std::vector<std::pair<SharedMember<Book>, double>> new_prof_books, extrap_books;

    for (auto &bc : library_on_market_) {
        if (bc.first->hasMarket()) {
            // The book is still on the market, so we'll have to extrapolate using profit stream
            // beliefs
            extrap_books.push_back(std::make_pair(bc.first, bc.second.get().quality));
        }
        else {
            // The book just left the market
            new_prof_books.push_back(std::make_pair(bc.first, bc.second.get().quality));
        }
    }

    // Also need to go through books that we don't own, using predicted quality
    for (auto &b : book_cache_market_) {
        if (b->hasMarket())
            extrap_books.push_back(std::make_pair(b, quality(b)));
        else
            new_prof_books.push_back(std::make_pair(b, quality(b)));
    }

    double weaken = creativity_->priorWeight();

    if (not new_prof_books.empty()) {
        MatrixXd X(new_prof_books.size(), profit_belief_->K());
        VectorXd y(new_prof_books.size());

        size_t i = 0;
        for (auto &bq : new_prof_books) {
            auto &book = bq.first;
            // Calculate the book's total profit
            y[i] = book->lifeRevenue() - book->marketPeriods() * cost_fixed - book->lifeSales() * cost_unit;
            X.row(i) = Profit::profitRow(bq.second, bq.first->order(), creativity_->market_books_lagged);
            i++;
        }

        *profit_belief_ = std::move(*profit_belief_).weaken(weaken).update(y, X);
    }
    else if (weaken != 1.0) {
        *profit_belief_ = std::move(*profit_belief_).weaken(weaken);
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
        for (auto &bq : extrap_books) {
            auto &book = bq.first;
            // Calculate the book's total profit (ignoring initial creation cost)
            y[i] = book->lifeRevenue() - book->marketPeriods() * cost_fixed - book->lifeSales() * cost_unit;
            X.row(i) = Profit::profitRow(bq.second, bq.first->order(), creativity_->market_books_lagged);
            i++;
        }

        // extrapolation uses just-updated non-extrapolation as prior, with no weakening
        profit_belief_extrap_.reset(new Profit(profit_belief_->update(y, X)));
    }
}

void Reader::intraInitialize() {
    auto nb = creativity_->newBooks();
    for (auto &bm : nb.first) {
        book_cache_.insert(bm);
        book_cache_market_.insert(bm);
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

        bool for_sale = book->hasMarket();
        bool shared = false;
        if (creativity_->sharing()) {
            for (auto &f : friends()) {
                if (f->library().count(book) > 0) {
                    // My friend owns the book
                    if (book->author() != f or not book->hasMarket()) {
                        // ... and either isn't the author, or is (and the book has left the market)
                        shared = true;
                        break;
                    }
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
            std::shuffle(bnu.second.begin(), bnu.second.end(), Random::rng());
    }

    lock.write(); // Time to get serious.

    for (auto &del : cache_del) book_cache_.erase(del);

    std::set<SharedMember<Book>> new_books;
    double money = assets()[creativity_->money];
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
                reserved_books_.emplace(book, false);
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
                reserved_books_.emplace(book, true);
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
        res->buy();
    }
    // Also remove any cost associated with obtaining pirated copies of books:
    assets()[creativity_->money] -= reserved_piracy_cost_;
    reserved_piracy_cost_ = 0.0;
    reservations_.clear();

    // Reset the "new" books list to the set of just-bought and just-pirated books
    library_new_.clear();

    for (auto &new_book : reserved_books_) {
        auto &book = new_book.first;
        auto &pirated = new_book.second;

        auto inserted = library_.emplace(
                SharedMember<Book>(book),
                BookCopy(book->qualityDraw(), pirated ? BookCopy::Status::pirated : BookCopy::Status::purchased, simulation()->t()));

        // If the book is still on the market (which it must be if we just bought it, and might be
        // if we just pirated it), stash a copy in library_no_market_ so that it will (eventually)
        // get incorporated into profit stream beliefs.
        if (book->hasMarket())
            library_on_market_.emplace(book, std::ref(inserted.first->second));

        // And it's always considered a "new" (to us) book, so stash it there:
        library_new_.emplace(book, std::ref(inserted.first->second));

        // If we obtained a pirated book, record that in the base book metadata:
        if (pirated) book->recordPiracy(1);
        //else {} -- // book->recordSale() not needed here: it's called during the BookMarket transaction
    }
    reserved_books_.clear();

    // "Eat" any money leftover
    double money = assets().remove(creativity_->money);
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
