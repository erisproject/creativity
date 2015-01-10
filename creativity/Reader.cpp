#include "creativity/Creativity.hpp"
#include "creativity/Reader.hpp"
#include "creativity/Book.hpp"
#include "creativity/BookMarket.hpp"
#include <eris/Random.hpp>
#include <algorithm>
#include <map>
#include <random>
#include <cmath>

using namespace eris;
using namespace Eigen;

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
    profit_belief_(creativity_->parameters.dimensions),
    profit_belief_extrap_(profit_belief_),
    demand_belief_(creativity_->parameters.dimensions),
    quality_belief_(belief::Quality::parameters()),
    profit_stream_beliefs_{{1, 1}}
{}

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
    // Update the various profit, demand, and quality beliefs
    updateBeliefs();

    double income_available = assets()[creativity_->money] + income;

    auto sim = simulation();
    const double market_books = sim->countMarkets<BookMarket>();
    const double previous_books = wrote().size();

    // Any books not in new_prices_ but that were on the market last period will be removed from the
    // market
    new_prices_.clear();

    auto &rng = Random::rng();

    if (usableBelief(demand_belief_)) {
        // Figure out the expected profitability of each book; positive ones will be kept on the market
        std::map<double, std::vector<std::pair<SharedMember<Book>, double>>, std::greater<double>> profitability;
        for (auto &book : wrote_market_) {
            auto book_mkt = book->market();

            auto max = demand_belief_.argmaxP(book->quality(), book->lifeSales(), previous_books, market_books, cost_unit);
            const double &p = max.first;
            const double &q = max.second;
            const double profit = (p - cost_unit) * q - cost_fixed;

            if (profit > 0) {
                // Profitable to keep this book on the market, so do so
                profitability[profit].push_back(std::make_pair(book, p));
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
        // Bypass market predictions if our market demand belief is noninformative, or not
        // sufficiently informed, and use initial model action parameters instead.
        if (creativity_->parameters.initial.prob_keep > 0) {
            std::bernoulli_distribution keep(creativity_->parameters.initial.prob_keep);
            for (auto on_market = wrote_market_.cbegin(); income_available >= cost_fixed and on_market != wrote_market_.cend(); on_market++) {
                auto &book = *on_market;
                if (keep(rng)) {
                    income_available -= cost_fixed;
                    double new_price = (book->price() - cost_unit) * creativity_->parameters.initial.keep_price + cost_unit;
                    new_prices_.emplace(book, new_price);
                }
            }
        }
    }

    //// NEW BOOK CREATION:
    create_ = false;
    if (income_available >= cost_fixed) {
        // Creating a book requires an ability to predict profit (to determine whether creation is
        // worthwhile), and an ability to predict demand (to determine the initial price)
        if (usableBelief(profit_belief_) and usableBelief(demand_belief_)) { // NB: profit_belief_extrap_ is a copy or update of profit_belief_, so will also be usable
            // Create a book if E(profit) > effort required

            // Find the l that maximizes creation profit given current wealth (which would have come from
            // past sales, if non-zero) plus the fixed income we're about to receive, minus whatever we
            // decided to spend above to keep books on the market.

            // FIXME: is this right, w.r.t. Bayesian MC prediction?  (Is averaging being done too
            // early?)
            double effort = profit_belief_extrap_.argmaxL(
                    [this] (double l) -> double { return creationQuality(l); },
                    previous_books, market_books, income_available - cost_fixed
                    );
            double quality = creationQuality(effort);

            double exp_profit = profit_belief_extrap_.predict(quality, previous_books, market_books);

            // If the optimal effort level gives a book with positive expected profits, do it:
            if (exp_profit > effort) {
                // We're going to create, so calculate the optimal first-period price.  It's possible that
                // we get back cost_unit (and so predicted profit is non-positive); write the book anyway:
                // perhaps profits are expected to come in later periods?
                auto max = demand_belief_.argmaxP(quality, 0, previous_books, market_books, cost_unit);
                create_price_ = max.first;
                create_quality_ = quality;
                create_effort_ = effort;
                create_ = true;
            }
        }
        else {
            // If we have no useful profit belief yet, just use the initial values:
            if (creativity_->parameters.initial.prob_write > 0 and std::bernoulli_distribution(creativity_->parameters.initial.prob_write)(rng)) {
                double q = std::uniform_real_distribution<double>(creativity_->parameters.initial.q_min, creativity_->parameters.initial.q_max)(rng);
                double effort = creationEffort(q);
                // Make sure the required effort doesn't exceed the available funds
                if (income_available >= cost_fixed + effort) {
                    create_ = true;
                    create_price_ = cost_unit + std::uniform_real_distribution<double>(creativity_->parameters.initial.p_min, creativity_->parameters.initial.p_max)(rng);
                    create_quality_ = q;
                    create_effort_ = effort;
                }
            }
        }
    }
}

void Reader::interApply() {
    // Give potential income (this has to be before authorship decision, since authorship requires
    // giving up some potential income: actual income will be this amount minus whatever is given up
    // to author)
    assets()[creativity_->money] += income;

    auto sim = simulation();
    SharedMember<Book> newbook;
    if (create_) {
        // The cost (think of this as an opportunity cost) of creating, and the first period fixed
        // cost:
        assets()[creativity_->money] -= create_effort_ + cost_fixed;

        // The book is centered at the reader's position, plus some noise we add below
        Position bookPos{position()};

        auto qdraw = [this] (const Book &book) -> double {
            return std::max(0.0, book.quality() + writer_quality_sd * Random::rstdnorm());
        };
        newbook = sim->spawn<Book>(creativity_, bookPos, sharedSelf(), wrote_.size(), create_price_, create_quality_, qdraw);

        /// If enabled, add some noise in a random direction to the position
        if (writer_book_sd > 0) {
            double step_dist = writer_book_sd * Random::rstdnorm();
            newbook->moveBy(step_dist * Position::random(bookPos.dimensions));
        }

        wrote_.insert(wrote_.end(), newbook);
        wrote_market_.insert(newbook);
        auto ins = library_.emplace(SharedMember<Book>(newbook), BookCopy(create_quality_, BookCopy::Status::wrote, sim->t()));
        library_unlearned_.emplace(newbook, std::ref(ins.first->second));
    }

    std::vector<SharedMember<Book>> remove;
    for (auto &b : wrote_market_) {
        // If it's staying on the market, update the price and incur the fixed cost
        if (new_prices_.count(b) > 0) {
            b->market()->setPrice(new_prices_[b]);
            assets()[creativity_->money] -= cost_fixed;
        }
        else if (not create_ or b != newbook) {
            // No new price, which means we're removing the book from the market
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
        // FIXME: check for proper Bayesian MC averaging
        auto &q_b = const_cast<belief::Quality&>(quality_belief_);
        q_hat = q_b.predict(b);
    }
    else {
        // No informative beliefs above quality, so just use initial parameter mean (= midpoint)
        q_hat = (creativity_->parameters.initial.q_max + creativity_->parameters.initial.q_min) / 2.0;
    }
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

const belief::Profit& Reader::profitBelief() const { return profit_belief_; }
const belief::Profit& Reader::profitExtrapBelief() const { return profit_belief_extrap_; }
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

bool Reader::usableBelief(const belief::Linear &model) const {
    return not(model.noninformative()) and model.n() - model.K() >= creativity_->parameters.initial.belief_threshold;
}

void Reader::updateBeliefs() {
    updateQualityBelief();
    updateDemandBelief();
    updateProfitStreamBelief();
    updateProfitBelief();

    // Clear any "on-market" book references that aren't on the market anymore, since they just got
    // incorporated into the beliefs above.
    std::vector<SharedMember<Book>> remove;
    for (auto &book : library_unlearned_) {
        if (not book.first->hasMarket())
            remove.push_back(book.first);
    }
    for (auto &book : remove)
        library_unlearned_.erase(book);
}

void Reader::updateQualityBelief() {
    double weaken = creativity_->priorWeight();
    if (weaken != 1.0) quality_belief_ = std::move(quality_belief_).weaken(weaken);

    // If we obtained any new books, update the demand belief with them
    if (not newBooks().empty()) {
        std::vector<SharedMember<Book>> books;
        VectorXd y(newBooks().size());
        size_t i = 0;
        for (auto &a : newBooks()) {
            books.push_back(a.first);
            y[i++] = a.second.get().quality;
        }

        quality_belief_ = std::move(quality_belief_).update(y, quality_belief_.bookData(books));
    }
}
void Reader::updateDemandBelief() {
    double weaken = creativity_->priorWeight();
    if (weaken != 1.0) demand_belief_ = std::move(demand_belief_).weaken(weaken);

    if (not newBooks().empty()) {
        std::vector<SharedMember<Book>> books;
        VectorXd y(newBooks().size());
        MatrixXd X(newBooks().size(), demand_belief_.K());
        // NB: this runs in the interoptimizer, which means t has already been incremented
        auto last_t = simulation()->t() - 1;
        size_t i = 0;
        for (auto &b : newBooks()) {
            y[i] = b.first->sales(last_t);
            X.row(i) = demand_belief_.bookRow(b.first, b.second.get().quality);
            i++;
        }

        demand_belief_ = std::move(demand_belief_).update(y, X);
    }
}
void Reader::updateProfitStreamBelief() {
    // Map book age into lists of books that survived on the market at least that long.  Only books
    // that have just left the market are considered (whether or not this reader bought or obtained
    // by piracy).
    std::map<unsigned int, std::vector<SharedMember<Book>>> to_learn;
    for (auto &bu : library_unlearned_) {
        auto &book = bu.first;
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
            profit_stream_beliefs_.emplace(age, age);
        }

        // Update the belief (existing or just-created) with the new data
        profit_stream_beliefs_.at(age) = profit_stream_beliefs_.at(age).update(y, X);
    }
}

void Reader::updateProfitBelief() {
    std::vector<std::pair<SharedMember<Book>, std::reference_wrapper<BookCopy>>> new_prof_books, extrap_books;
    for (auto &bc : library_unlearned_) {
        if (bc.first->hasMarket()) {
            // The book is still on the market, so we'll have to extrapolate using profit stream
            // beliefs
            extrap_books.push_back(bc);
        }
        else {
            // The book either just left the market or we just obtained an off-market book via
            // piracy; in either case, incorporate it into the profit belief.
            new_prof_books.push_back(bc);
        }
    }

    double weaken = creativity_->priorWeight();
    if (weaken != 1.0) profit_belief_ = std::move(profit_belief_).weaken(weaken);

    if (not new_prof_books.empty()) {
        MatrixXd X(new_prof_books.size(), profit_belief_.K());
        VectorXd y(new_prof_books.size());

        size_t i = 0;
        for (auto &bc : new_prof_books) {
            auto &book = bc.first;
            // Calculate the book's total (optimal) profit, ignoring trailing negative profit
            // periods (which would be non-optimal for this reader (and probably for the actual
            // author, though not necessarily because the actual author might have different
            // costs)).
            double profit_total = 0;
            for (auto t = book->outOfPrint()-1; t >= book->created(); t--) {
                double prof_t = book->revenue(t) - cost_fixed - cost_unit * book->sales(t);
                profit_total += prof_t;
                if (profit_total < 0) {
                    // Total profit is negative, which means that keeping the book past this point
                    // was non-optimal, so reset to 0 (i.e. don't count profits from this point on).
                    // This only happens when the period's profits are negative, but doesn't
                    // *necessarily* happen in such a case: it could also be that negative profits
                    // are followed by larger, positive profits.
                    profit_total = 0;
                }
            }
            y[i] = profit_total;
            X.row(i) = profit_belief_.profitRow(bc.first, bc.second.get().quality);
            i++;
        }

        profit_belief_ = profit_belief_.update(y, X);
    }

    // Extrapolate based on profit stream predictions for profit levels for books that are still on
    // the market (since their profit level is not yet finalized)
    if (extrap_books.empty()) {
        // No extrapolation books, so the "extrapolation" belief is just the profit belief
        profit_belief_extrap_ = profit_belief_;
    }
    else {
        MatrixXd X(extrap_books.size(), profit_belief_.K());
        VectorXd y(extrap_books.size());

        size_t i = 0;
        for (auto &bc : extrap_books) {
            // Look for the largest model that doesn't exceed the book's age, then use it for
            // prediction
            for (auto it = profit_stream_beliefs_.rbegin(); it != profit_stream_beliefs_.rend(); it++) {
                if (it->first <= bc.first->age() and it->second.n() >= 1) {
                    // We have a winner:
                    y[i] = it->second.predict(bc.first);
                    X.row(i) = profit_belief_.profitRow(bc.first, bc.second.get().quality);
                    i++;
                    break;
                }
            }
        }

        // NB: extrapolation uses just-updated non-extrapolation as prior, with no weakening
        profit_belief_extrap_ = profit_belief_.update(y, X);
    }
}

void Reader::intraInitialize() {
    auto nb = creativity_->newBooks();
    for (auto &bm : nb.first) {
        book_cache_.insert(bm);
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
                    if (book->author() != f) {
                        // ... and isn't the author (who presumably wouldn't share it)
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
        // if we just pirated it), stash a copy in library_unlearned_ so that it will (eventually)
        // get incorporated into profit stream beliefs.
        if (book->hasMarket())
            library_unlearned_.emplace(book, std::ref(inserted.first->second));

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
