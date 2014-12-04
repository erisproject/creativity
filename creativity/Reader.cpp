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
        belief::Demand &&demand, belief::Profit &&profit, belief::Quality &&quality,
        double cFixed, double cUnit, double inc
        )
    : WrappedPositional<agent::AssetAgent>(pos, {-creativity->parameters.boundary, -creativity->parameters.boundary}, {creativity->parameters.boundary, creativity->parameters.boundary}),
    cost_fixed{cFixed}, cost_unit{cUnit}, income{inc},
    creativity_{std::move(creativity)},
    profit_belief_{std::move(profit)}, profit_belief_extrap_{profit_belief_},
    demand_belief_{std::move(demand)}, quality_belief_{std::move(quality)},
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

const double& Reader::u() const {
    return u_curr_;
}
const double& Reader::uLifetime() const {
    return u_lifetime_;
}
const std::unordered_set<SharedMember<Book>>& Reader::newBooks() const {
    return library_new_;
}
const std::unordered_set<SharedMember<Book>>& Reader::newPurchased() const {
    return library_new_purchased_;
}
const std::unordered_set<SharedMember<Book>>& Reader::newPirated() const {
    return library_new_pirated_;
}
const std::unordered_set<SharedMember<Book>>& Reader::libraryPirated() const {
    return library_pirated_;
}
const std::unordered_set<SharedMember<Book>>& Reader::libraryPurchased() const {
    return library_purchased_;
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
    if (effort < 0)
        throw std::domain_error("Reader::creationQuality() error: effort cannot be negative");
    if (creation_shape < 0)
        throw std::logic_error("Reader::creationQuality() error: creation_shape cannot be negative");
    if (creation_scale < 0)
        throw std::logic_error("Reader::creationQuality() error: creation_scale cannot be negative");
    if (creation_scale == 0) return 0.0;

    return creation_scale * (
            creation_shape == 1.0 ? std::log(effort + 1.0) : // Special case for log
            creation_shape == 0.5 ? 2.0 * (std::sqrt(effort + 1.0) - 1.0) : // Slightly more efficient sqrt calculation
            creation_shape == 0.0 ? effort : // Simple linear
            (std::pow(effort + 1.0, 1.0 - creation_shape) - 1.0) / (1 - creation_shape) // Otherwise use the full formula
    );
}

double Reader::creationEffort(double quality) const {
    if (quality < 0)
        throw std::domain_error("Reader::creationEffort() error: quality cannot be negative");
    if (creation_shape < 0)
        throw std::logic_error("Reader::creationEffort() error: creation_shape cannot be negative");
    if (creation_scale < 0)
        throw std::logic_error("Reader::creationEffort() error: creation_scale cannot be negative");
    if (creation_scale == 0) return 0.0;

    return
        creation_shape == 1.0 ? std::exp(quality / creation_scale) - 1 : // The special log case
        creation_shape == 0.0 ? quality / creation_scale : // Simple linear case
        std::pow(1.0 + quality * (1.0-creation_shape) / creation_scale, 1.0/(1.0 - creation_shape)) - 1;
}

double Reader::piracyCost() const {
    return creativity_->parameters.cost_unit;
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

    // Bypass market predictions if our market demand belief is noninformative:
    if (demand_belief_.noninformative()) {
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
    else {
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

    //// NEW BOOK CREATION:
    create_ = false;
    if (income_available >= cost_fixed) {
        if (profit_belief_.noninformative()) {
            // If we have no useful profit belief yet, just use the initial values:
            if (creativity_->parameters.initial.prob_write > 0 and std::bernoulli_distribution(creativity_->parameters.initial.prob_write)(rng)) {
                double q = std::uniform_real_distribution<double>(creativity_->parameters.initial.q_min, creativity_->parameters.initial.q_max)(rng);
                double effort = creationEffort(q);
                // Make sure the required effort doesn't exceed the available funds
                if (income_available >= cost_fixed + effort) {
                    create_ = true;
                    create_price_ = std::uniform_real_distribution<double>(creativity_->parameters.initial.p_min, creativity_->parameters.initial.p_max)(rng);
                    create_quality_ = q;
                    create_effort_ = effort;
                    ERIS_DBG(id() << " creating (via initial parameters) at q=" << q << ", l=" << effort <<", p=" << create_price_);
                }
            }
        }
        else {
            // Otherwise, create a book if E(profit) > effort required

            // Find the l that maximizes creation profit given current wealth (which would have come from
            // past sales, if non-zero) plus the fixed income we're about to receive, minus whatever we
            // decided to spend above to keep books on the market.
            double l_max = profit_belief_extrap_.argmaxL(
                    [this] (double l) -> double { return creationQuality(l); },
                    previous_books, market_books, assets()[creativity_->money] + income - cost_fixed
                    );
            ERIS_DBGVAR(id());
            ERIS_DBGVAR(l_max);
            double quality = creationQuality(l_max);
            ERIS_DBGVAR(quality);

            double exp_profit = profit_belief_.predict(quality, previous_books, market_books);
            ERIS_DBGVAR(exp_profit);

            // If the optimal effort level gives a book with positive expected profits, do it:
            if (l_max < exp_profit) {
                ERIS_DBG("CREATING! =)");
                // We're going to create, so calculate the optimal first-period price.  It's possible that
                // we get back cost_unit (and so predicted profit is non-positive); write the book anyway:
                // perhaps profits are expected to come in later periods?
                auto max = demand_belief_.argmaxP(quality, 0, previous_books, market_books, cost_unit);
                create_price_ = max.first;
                create_quality_ = quality;
                create_effort_ = l_max;
                create_ = true;
            }
            else ERIS_DBG("not creating. :(");
        }
    }
}

void Reader::interApply() {
    // Give potential income (this has to be before authorship decision, since authorship requires
    // giving up some potential income: actual income will be this amount minus whatever is given up
    // to author)
    assets()[creativity_->money] += income;

    SharedMember<Book> newbook;
    if (create_) {
        // The cost (think of this as an opportunity cost) of creating:
        assets()[creativity_->money] -= create_effort_ + cost_fixed;

        // The book is centered at the reader's position, plus some noise we add below
        Position bookPos{position()};

        auto qdraw = [this] (const Book &book, const Reader &) -> double {
            return std::max(0.0, book.quality() + writer_quality_sd * stdnormal(Random::rng()));
        };
        newbook = simulation()->spawn<Book>(creativity_, bookPos, sharedSelf(), wrote_.size(), create_price_, create_quality_, qdraw);

        /// If enabled, add some noise in a random direction to the position
        if (writer_book_sd > 0) {
            double step_dist = std::normal_distribution<double>(0, writer_book_sd)(Random::rng());
            newbook->moveBy(step_dist * Position::random(bookPos.dimensions));
        }

        wrote_.insert(wrote_.end(), newbook);
        wrote_market_.insert(newbook);
        library_.emplace(newbook, create_quality_);
        library_unlearned_.insert(newbook);
    }

    std::vector<SharedMember<Book>> remove;
    for (auto &b : wrote_market_) {
        if (new_prices_.count(b) > 0) {
            // Set the new price
            b->market()->setPrice(new_prices_[b]);
            assets()[creativity_->money] -= cost_fixed;
        }
        else if (not create_ or b != newbook) {
            // No new price, which means we remove the book from the market
            remove.push_back(b);
        }
    }

    for (auto &b : remove) {
        wrote_market_.erase(b);
        simulation()->remove(b->market());
    }

    // Finally, move N(0,0.25) in a random direction
    double step_dist = std::normal_distribution<double>(0, 0.25)(Random::rng());
    moveBy(step_dist * Position::random(position().dimensions));
}

double Reader::uBook(SharedMember<Book> b) const {
    double u = evalPolynomial(distance(b), u_poly_);
    u += quality(b);
    if (u < 0) u = 0.0;
    return u;
}

double Reader::quality(SharedMember<Book> b) const {
    auto found = library_.find(b);
    if (found != library_.end())
        return found->second;

    found = quality_predictions_.find(b);
    if (found != quality_predictions_.end())
        return found->second;

    double q_hat;
    if (quality_belief_.noninformative()) {
        // No informative beliefs above quality, so just use initial parameter midpoint
        q_hat = (creativity_->parameters.initial.q_max + creativity_->parameters.initial.q_min) / 2.0;
    }
    else {
        // Otherwise we need to use the quality belief to predict the quality
        auto &q_b = const_cast<belief::Quality&>(quality_belief_);
        q_hat = q_b.predict(b);
    }
    quality_predictions_.emplace(b, q_hat);

    return q_hat;
}

double Reader::penalty(unsigned long n) const {
    return evalPolynomial(n, pen_poly_);
}

const std::unordered_map<SharedMember<Book>, double>& Reader::library() const { return library_; }

void Reader::receiveProfits(SharedMember<Book> book, const Bundle &revenue) {
    assets() += revenue;
    assets()[creativity_->money] -= book->currSales() * cost_unit;
}

const belief::Profit& Reader::profitBelief() const { return profit_belief_; }
const belief::Profit& Reader::profitExtrapBelief() const { return profit_belief_extrap_; }
const belief::Demand& Reader::demandBelief() const { return demand_belief_; }
const belief::Quality& Reader::qualityBelief() const { return quality_belief_; }
const belief::ProfitStream& Reader::profitStreamBelief(unsigned int age) const {
    // Get the first belief > age
    auto it = profit_stream_beliefs_.upper_bound(age);
    if (it == profit_stream_beliefs_.begin())
        // This probably means age=0 got passed in, since profit_stream_beliefs_ starts with a 1.
        throw std::runtime_error("Invalid age (" + std::to_string(age) + ") passed to Reader::profitStreamBelief");
    // upper_bound gives first greater than the given value, so we need to back up one:
    it--;
    return it->second;
}
const std::map<unsigned int, belief::ProfitStream>& Reader::profitStreamBeliefs() const {
    return profit_stream_beliefs_;
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
        if (not book->hasMarket())
            remove.push_back(book);
    }
    for (auto &book : remove)
        library_unlearned_.erase(book);
}

void Reader::updateQualityBelief() {
    // If we obtained any new books, update the demand belief with them
    if (not newBooks().empty()) {
        std::vector<SharedMember<Book>> books;
        VectorXd y(newBooks().size());
        size_t i = 0;
        for (auto &a : newBooks()) {
            books.push_back(a);
            y[i++] = library_[a];
        }

        quality_belief_ = quality_belief_.update(y, quality_belief_.bookData(books));
    }
}
void Reader::updateDemandBelief() {
    if (not newBooks().empty()) {
        std::vector<SharedMember<Book>> books;
        VectorXd y(newBooks().size());
        MatrixXd X(newBooks().size(), demand_belief_.K());
        auto t = simulation()->t();
        size_t i = 0;
        for (auto &b : newBooks()) {
            y[i] = b->sales(t);
            X.row(i) = demand_belief_.bookRow(b, library_[b]);
            i++;
        }

        demand_belief_ = demand_belief_.update(y, X);
    }
}
void Reader::updateProfitStreamBelief() {
    // Map book age into lists of books that survived on the market at least that long.  Only books
    // that have left the market are considered.  Note that being added here doesn't mean the book
    // *just* left the market: it might have just left, but it also might be a pirated book that
    // left the market a long time ago.
    std::map<unsigned int, std::vector<SharedMember<Book>>> to_learn;
    for (auto &book : library_unlearned_) {
        if (not book->hasMarket()) {
            unsigned long periods = book->marketPeriods();
            // We only look for the predefined age values, and only include books with
            // marketPeriods() strictly greater than the predefined age values (because a book on
            // the market for x periods can only contribute to models with (x-1) past profit
            // variables.
            //
            // However, it's possible that the book was kept on the market for too long (i.e. it had
            // no sales on the last 2+ periods), so reduce `periods` until either the last or
            // second-last referenced period had positive sales.  (The last having zero sales is
            // fine, because it will be the dependent variable for the model including the
            // second-last).
            const unsigned long &created = book->created();
            while (periods > 1 and book->sales(created + periods - 1) == 0 and book->sales(created + periods - 2) == 0)
                periods--;

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

    for (auto &learn : to_learn) {
        const auto &age = learn.first;
        RowVectorXd y(learn.second.size());
        MatrixXd X(learn.second.size(), age);

        size_t row = 0;
        for (auto &book : learn.second) {
            double cumul_rev = 0;
            const unsigned long &created = book->created();
            for (unsigned int i = 0; i < age; i++) {
                double r_i = book->revenue(created + i);
                X(row, i) = r_i;
                cumul_rev += r_i;
            }
            y[row] = book->lifeRevenue() - cumul_rev;
            row++;
        }

        // Create a new belief of the given size if none exists yet
        if (profit_stream_beliefs_.count(age) == 0) {
            // We'll just use a highly noninformative prior, so the update below will determine
            // almost everything.
            VectorXd beta = VectorXd::Zero(age, 1);

            profit_stream_beliefs_.emplace(age, age);
        }

        // Update the belief (existing or just-created) with the new data
        profit_stream_beliefs_.at(age) = profit_stream_beliefs_.at(age).update(y, X);
    }
}

void Reader::updateProfitBelief() {
    std::vector<SharedMember<Book>> new_prof_books, extrap_books;
    for (auto &book : library_unlearned_) {
        if (book->hasMarket()) {
            // The book is still on the market, so we'll have to extrapolate using profit stream
            // beliefs
            extrap_books.push_back(book);
        }
        else {
            // The book either just left the market or we just obtained an off-market book via
            // piracy; in either case, incorporate it into the profit belief.
            new_prof_books.push_back(book);
        }
    }

    if (not new_prof_books.empty()) {
        MatrixXd X(new_prof_books.size(), profit_belief_.K());
        VectorXd y(new_prof_books.size());

        size_t i = 0;
        for (auto &book : new_prof_books) {
            y[i] = book->lifeRevenue();
            X.row(i) = profit_belief_.profitRow(book, library_[book]);
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
        for (auto &book : extrap_books) {
            // Look for the largest model that doesn't exceed the book's age, then use it for
            // prediction
            for (auto it = profit_stream_beliefs_.rbegin(); it != profit_stream_beliefs_.rend(); it++) {
                if (it->first <= book->age() and it->second.n() >= 1) {
                    // We have a winner:
                    y[i] = it->second.predict(book);
                    X.row(i) = profit_belief_.profitRow(book, library_[book]);
                    i++;
                    break;
                }
            }
        }

        // NB: extrapolation uses just-updated non-extrapolation as prior
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
    for (auto &res : reservations_) {
        res->buy();
    }
    assets()[creativity_->money] -= reserved_piracy_cost_;
    reserved_piracy_cost_ = 0.0;
    reservations_.clear();

    library_new_.clear();
    library_new_pirated_.clear();
    library_new_purchased_.clear();

    for (auto &new_book : reserved_books_) {
        auto &book = new_book.first;
        auto &pirated = new_book.second;

        double q = book->qualityDraw(*this);
        library_.emplace(book, q);
        library_unlearned_.insert(book);
        library_new_.insert(book);
        if (pirated) {
            book->recordPiracy(1);
            library_pirated_.insert(book);
            library_new_pirated_.insert(book);
        }
        else {
            // book->recordSale not needed here: it's called during the BookMarket transaction
            library_purchased_.insert(book);
            library_new_purchased_.insert(book);
        }
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
