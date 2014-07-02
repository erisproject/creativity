#include "creativity/Reader.hpp"
#include "creativity/Book.hpp"
#include "creativity/BookMarket.hpp"
#include <eris/Random.hpp>
#include <algorithm>
#include <map>

using namespace eris;
using namespace Eigen;

namespace creativity {

const std::vector<double> Reader::default_polynomial{{4., -1.}};
const std::vector<double> Reader::default_penalty_polynomial{{0, 0, 0.25}};
const std::vector<unsigned long> Reader::profit_stream_ages_{{1,2,3,4,5,7,10}};

Reader::Reader(
        const Position &pos, const Position &b1, const Position &b2,
        belief::Demand &&demand, belief::Profit &&profit, belief::Quality &&quality,
        const double &cFixed, const double &cUnit, const double &income
        )
    : WrappedPositional<agent::AssetAgent>(pos, b1, b2),
    profit_belief_{std::move(profit)}, profit_belief_extrap_{profit_belief_},
    demand_belief_{std::move(demand)}, quality_belief_{std::move(quality)},
    c_fixed_{cFixed}, c_unit_{cUnit}, income_{income}
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

const std::vector<SharedMember<Book>>& Reader::wrote() const {
    return wrote_;
}

const std::vector<double>& Reader::uPolynomial() const {
    return u_poly_;
}

double Reader::evalPolynomial(const double &x, const std::vector<double> &polynomial) {
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

double Reader::creationQuality(const double &effort) const {
    return creation_coefs[0] + creation_coefs[1] * std::pow(effort, creation_coefs[2]);
}

void Reader::interOptimize() {
    // Update the various profit, demand, and quality beliefs
    updateBeliefs();

    // Clear any "on-market" book references that aren't on the market anymore
    auto it = library_market_.begin();
    while (it != library_market_.end()) {
        if (not (*it)->hasMarket())
            it = library_market_.erase(it);
        else
            it++;
    }

    double income_available = assets()[MONEY] + income_;

    auto sim = simulation();
    const double market_books = sim->countMarkets<BookMarket>();
    const double previous_books = wrote().size();

    // Figure out the expected profitability of each book; positive ones will be kept on the market
    std::map<double, std::vector<std::pair<SharedMember<Book>, double>>, std::greater<double>> profitability;
    for (auto &book : wrote_market_) {
        auto book_mkt = book->market();

        auto max = demand_belief_.argmaxP(book->quality(), book->lifeSales(), previous_books-1, market_books, c_unit_);
        const double &p = max.first;
        const double &q = max.second;
        const double profit = (p - c_unit_) * q - c_fixed_;

        if (profit > 0) {
            // Profitable to keep this book on the market, so do so
            profitability[profit].push_back(std::make_pair(book, p));
        }
    }

    // Store the decided prices for each book staying on the market in decreasing order of
    // profitability UNLESS we can't afford the fixed cost of keeping the book on the market.
    // Anything with a negative expected profit or costs that can't be covered will be removed from
    // the market.  Note that the costs calculated here aren't actually incurred until interApply().
    new_prices_.clear();
    while (income_available >= c_fixed_ and not profitability.empty()) {
        auto &books = profitability.begin()->second;
        if (income_available < c_fixed_ * books.size())
            // We don't have enough to do all the books at this profitability level, so shuffle them
            // so that we choose randomly
            std::shuffle(books.begin(), books.end(), Random::rng());

        for (auto &b : books) {
            if (income_available < c_fixed_) break;
            income_available -= c_fixed_;
            new_prices_.emplace(b.first, b.second);
        }
    }
    // Any books not in new_prices_ but that were on the market last period will be removed from the
    // market


    //// NEW BOOK CREATION:
    // Create a book if E(profit) > effort required
    create_ = false;

    // Find the l that maximizes creation profit given current wealth (which would have come from
    // past sales, if non-zero) plus the fixed income we're about to receive, minus whatever we
    // decided to spend above to keep books on the market.
    if (income_available >= c_fixed_) {
        double l_max = profit_belief_extrap_.argmaxL(
                [this] (const double &l) -> double { return creationQuality(l); },
                previous_books, market_books, assets()[MONEY] + income_ - c_fixed_
                );
        double quality = creationQuality(l_max);

        double exp_profit = profit_belief_.predict(quality, previous_books, market_books);

        // If the optimal effort level gives a book with positive expected profits, do it:
        if (l_max < exp_profit) {
            // We're going to create, so calculate the optimal first-period price.  It's possible that
            // we get back c_unit_ (and so predicted profit is non-positive); write the book anyway:
            // perhaps profits are expected to come in later periods?
            auto max = demand_belief_.argmaxP(quality, 0, wrote_.size()-1, market_books, c_unit_);
            create_price_ = max.first;
            create_quality_ = quality;
            create_effort_ = l_max;
            create_ = true;
        }
    }
}

void Reader::interApply() {
    // Give potential income (this has to be before authorship decision, since authorship requires
    // giving up some potential income: actual income will be this amount minus whatever is given up
    // to author)
    assets() += {MONEY, income_};

    SharedMember<Book> newbook;
    if (create_) {
        // The cost (think of this as an opportunity cost) of creating:
        assets() -= {MONEY, create_effort_ + c_fixed_};

        // The book is centered at the reader's position, plus some noise we add below
        Position bookPos{position()};

        auto qdraw = [this] (const Book &book, const Reader &) -> double {
            return std::max(0.0, book.quality() + writer_quality_sd * stdnormal(Random::rng()));
        };
        newbook = simulation()->create<Book>(bookPos, sharedSelf(), wrote_.size(), create_price_, create_quality_, qdraw);

        /// If enabled, add some noise in a random direction to the position
        if (writer_book_sd > 0) {
            std::normal_distribution<double> step_dist{0, writer_book_sd};
            newbook->moveBy(step_dist(Random::rng()) * Position::random(bookPos.dimensions));
        }

        wrote_.push_back(newbook);
        wrote_market_.insert(newbook);
        library_.emplace(newbook, create_quality_);
        library_market_.insert(newbook);
    }

    std::vector<SharedMember<Book>> remove;
    for (auto &b : wrote_market_) {
        if (new_prices_.count(b) > 0) {
            // Set the new price
            b->market()->setPrice(new_prices_[b]);
            assets() -= {MONEY, c_fixed_};
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
    std::normal_distribution<double> step_dist{0, 0.25};
    moveBy(step_dist(Random::rng()) * Position::random(position().dimensions));
}

double Reader::uBook(SharedMember<Book> b) const {
    double u = evalPolynomial(distance(b), u_poly_);
    double q_DEBUG = quality(b);
    //std::cerr << "Evaluating book " << b << ": quality=" << q_DEBUG << ", totalu=" << q_DEBUG + u << "\n";
    u += q_DEBUG;//quality(b);
    if (u < 0) u = 0.0;
    return u;
}

double Reader::quality(SharedMember<Book> b) const {
    auto found = library_.find(b);
    if (found != library_.end())
        return found->second;

    // Otherwise we need to use the quality belief to predict the quality
    return quality_belief_.predict(b);
}

double Reader::penalty(unsigned long n) const {
    return evalPolynomial(n, pen_poly_);
}

const std::unordered_map<SharedMember<Book>, double>& Reader::library() { return library_; }

void Reader::receiveProfits(SharedMember<Book> book, const Bundle &revenue) {
    assets() += revenue - Bundle{MONEY, book->currSales() * c_unit_};
}

const belief::Profit& Reader::profitBelief() { return profit_belief_; }
const belief::Profit& Reader::profitExtrapBelief() { return profit_belief_extrap_; }
const belief::Demand& Reader::demandBelief() { return demand_belief_; }
const belief::Quality& Reader::qualityBelief() { return quality_belief_; }
const belief::ProfitStream& Reader::profitStreamBelief(unsigned long age) {
    if (profit_stream_beliefs_.empty()) {
        profit_stream_beliefs_.emplace(std::piecewise_construct,
                std::tuple<unsigned long>(age),
                std::tuple<size_t>(age));
    }

    auto it = profit_stream_beliefs_.upper_bound(age);
    if (it == profit_stream_beliefs_.begin())
        throw std::runtime_error("Invalid age (" + std::to_string(age) + ") passed to Reader::profitStreamBelief");
    it--;
    return it->second;
}

void Reader::updateBeliefs() {
    updateQualityBelief();
    updateDemandBelief();
    updateProfitStreamBelief();
    updateProfitBelief();
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
    // Map age into lists of just-left-market books of at least that age
    std::map<unsigned long, std::vector<SharedMember<Book>>> just_left;
    ERIS_DBGVAR(library_market_.size());
    for (auto &book : library_market_) {
        if (not book->hasMarket()) {
            const unsigned long periods = book->marketPeriods();
            // We only look for the predefined age values, and only include books with
            // marketPeriods() strictly greater than the predefined age values (because a book on
            // the market for x periods can only contribute to models with (x-1) past profit
            // variables.
            //
            // TODO: this could potentially change to create new profit stream age models as we
            // encounter books of longer ages.  One tricky thing to consider: if we have a book of
            // age 6, should we really use, say, a n=1, age=6 model when we have, say, n=10, age=5
            // and n=100, age=4 models also available?  Choose one?  Weighted average?
            for (unsigned long a : profit_stream_ages_) {
                if (periods > a)
                    just_left[a].push_back(book);
                else
                    break;
            }
        }
    }

    ERIS_DBG("Updating profit stream beliefs with:");
    for (auto jl : just_left) {
        ERIS_DBG("    " << jl.first << ": " << jl.second.size() << " books");
    }
    ERIS_DBG(".");

    for (auto &jl : just_left) {
        const unsigned long &age = jl.first;
        RowVectorXd y(jl.second.size());
        MatrixXd X(jl.second.size(), age);

        size_t row = 0;
        for (auto &book : jl.second) {
            double cumul_rev = 0;
            unsigned long created = book->created();
            for (unsigned long i = 0; i < age; i++) {
                double r_i = book->revenue(created + i);
                X(row, i) = r_i;
                cumul_rev += r_i;
            }
            y[row] = book->lifeRevenue() - cumul_rev;
        }

        // Create a new belief of the given size if none exists yet
        if (profit_stream_beliefs_.count(age) == 0) {
            // We'll just make up a value that remaining profits are equal to the last period's
            // profits, with some made up variance values.  This is certainly wrong, but it'll also
            // be constructed with a weak prior which should quickly wash out from new observations.
            VectorXd beta = VectorXd::Zero(age, 1);
            double s2 = 1.0;
            MatrixXd V = MatrixXd::Identity(beta.rows(), beta.rows());
            double n = 1e-6; // So that this is very weak when used as a prior

            profit_stream_beliefs_.emplace(std::piecewise_construct,
                    std::tuple<unsigned long>(age),
                    std::tuple<size_t>(age));
        }

        // Update the belief (existing or just-created) with the new data
        profit_stream_beliefs_.at(age) = profit_stream_beliefs_.at(age).update(y, X);
    }
}

void Reader::updateProfitBelief() {
    std::vector<SharedMember<Book>> new_prof_books, extrap_books;
    for (auto &book : library_market_) {
        if (book->hasMarket()) {
            // The book is still on the market, so we'll have to extrapolate using profit stream
            // beliefs
            extrap_books.push_back(book);
        }
        else {
            // The book just left the market (it's still in library_market_, from the previous
            // period), so incorporate it into the profit belief.
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
    for (auto &bm : NEW_BOOKS) {
        book_cache_.insert(bm);
    }
}

void Reader::intraOptimize() {
    auto lock = readLock();

    std::vector<SharedMember<Book>> cache_del;
    // Map utility-minus-price values to sets of market ids so that we can pick the book(s) where
    // net utility (i.e. u-p) is highest.  NB: this is sorted by utility in highest-to-lowest order
    // for the purposes of iteration.
    std::map<double, std::vector<SharedMember<Book>>, std::greater<double>> book_net_u;
    for (auto &book : book_cache_) {
        if (not book->hasMarket()) { // Market removed: remove from cache and continue
            cache_del.push_back(book);
            continue;
        }
        if (library_.count(book) > 0) {
            // Already have the book: remove from cache and continue
            cache_del.push_back(book);
            continue;
        }
        double u = uBook(book);
        double p = book->market()->price();
        if (u >= p) {
            // Good: the book is available and its utility (before any penalty) exceeds the purchase price
            // Store the net utility of this book:
            book_net_u[u - p].push_back(book);
        }
    }

    // Shuffle the order of any books with the same net utilities so that we'll be choosing randomly
    // from equally-preferred options.
    for (auto &bnu : book_net_u) {
        if (bnu.second.size() > 1)
            std::shuffle(bnu.second.begin(), bnu.second.end(), Random::rng());
    }

    lock.write(); // Time to get serious.

    for (auto &del : cache_del) book_cache_.erase(del);

    std::set<SharedMember<Book>> buy;
    double money = assets()[MONEY];
    double u_curr = u(money, buy); // the "no-books" utility
    // Keep looking at books as long as the net utility from buying the next best book exceeds the
    // penalty that will be incurred.
    while (not book_net_u.empty() and money > 0) {
        // Pull off the best remaining book:
        auto &s = book_net_u.begin()->second;
        auto book = std::move(s.back());
        s.pop_back();
        // If we just pulled off the last book at the current net utility level, delete the level.
        if (s.empty()) book_net_u.erase(book_net_u.begin());

        auto bm = book->market();
        lock.add(bm);
        // Perform some safety checks:
        // - if we've already added the book to the set of books to buy, don't add it again. (This
        //   should only be possible if there are multiple BookMarkets selling the same book).
        // - if the BookMarket says the book is no longer feasible, ignore it. (This shouldn't
        //   happen, but just in case).
        // - if we can't afford it, skip it.
        auto pinfo = bm->price(1);
        if (buy.count(book) == 0 and pinfo.feasible and pinfo.total <= money) {
            buy.insert(book);
            double u_with_book = u(money - pinfo.total, buy);
            if (u_with_book < u_curr) {
                // The best book lowered utility, so we're done.
                buy.erase(book);
                lock.remove(bm);
                break;
            }
            // Otherwise the purchase is a good one, so make a reservation.
            reservations_.push_front(bm->reserve(sharedSelf(), 1, pinfo.total));
            reserved_books_.insert(book);
            u_curr = u_with_book;
            money -= pinfo.total;
        }
        lock.remove(bm);
    }
}
void Reader::intraApply() {
    auto lock = writeLock();
    for (auto &res : reservations_) {
        res->buy();
    }
    for (auto &book : reserved_books_) {
        double q = book->qualityDraw(*this);
        library_.emplace(book, q);
        library_market_.insert(book);
    }
    reservations_.clear();
    library_new_.clear();
    std::swap(library_new_, reserved_books_); // Clear away into new_books_

    // "Eat" any money leftover
    double money = assets().remove(MONEY);
    // Finalize utility
    u_curr_ = u(money, library_new_);
    u_lifetime_ += u_curr_;
}
void Reader::intraReset() {
    auto lock = writeLock();
    reservations_.clear();
    reserved_books_.clear();
}

}
