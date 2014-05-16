#include "creativity/Reader.hpp"
#include "creativity/Book.hpp"
#include "creativity/BookMarket.hpp"
#include <eris/Random.hpp>
#include <algorithm>
#include <map>

using namespace eris;

namespace creativity {

const std::vector<double> Reader::default_polynomial{{4., -1.}};
const std::vector<double> Reader::default_penalty_polynomial{{0, 0, 0.25}};

Reader::Reader(
        const Position &pos, const Position &b1, const Position &b2,
        belief::Demand &&demand, belief::Profit &&profit, belief::ProfitStream &&stream, belief::Quality &&quality,
        const double &cFixed, const double &cUnit, const double &income
        )
    : WrappedPositional<agent::AssetAgent>(pos, b1, b2),
    profit_belief_{std::move(profit)}, profit_belief_extrap_{profit_belief_}, prof_stream_belief_{std::move(stream)},
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
    return new_books_;
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

    auto sim = simulation();

    const double previous_books = wrote().size();
    const double market_books = sim->countMarkets<BookMarket>();

    // Set the price of all previous books for the upcoming period.  Remove from market if expected
    // profit is negative.
    new_prices_.clear();
    for (auto &book : wrote_on_market_) {
        auto book_mkt = book->market();

        auto max = demand_belief_.argmaxP(book->quality(), book->lifeSales(), previous_books-1, market_books, c_unit_);
        const double &p = max.first;
        const double &q = max.second;
        const double profit = (p - c_unit_) * q - c_fixed_;

        if (profit > 0) {
            // Profitable to keep this book on the market, so do so
            new_prices_.emplace(book, p);
        }
        // Otherwise the book will be removed from the market
    }

    // Create a book if E(profit) > effort required

    // Find the l that maximizes creation profit
    double l_max = profit_belief_extrap_.argmaxL(
            [this] (const double &l) -> double { return creationQuality(l); },
            previous_books, market_books, assets()[MONEY]
            );
    ERIS_DBGVAR(l_max);
    double quality = creationQuality(l_max);

    ERIS_DBGVAR(quality);

    if (quality >= 0) ERIS_DBGVAR(profit_belief_.predict(quality, previous_books, market_books));

    // If the optimal effort level gives a book with positive expected profits, do it:
    if (quality >= 0 and l_max < profit_belief_.predict(quality, previous_books, market_books)) {
        // We're going to create, so calculate the optimal first-period price.  It's possible that
        // we get back c_unit_ (and so predicted profit is non-positive); write the book anyway:
        // perhaps profits are expected to come in later periods.
        auto max = demand_belief_.argmaxP(quality, 0, wrote_.size()-1, market_books, c_unit_);
        create_price_ = max.first;
        create_quality_ = quality;
        create_effort_ = l_max;
        create_ = true;
    }
    else {
        create_ = false;
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
        assets() -= {MONEY, create_effort_};

        // The book is centered at the reader's position:
        Position bookPos{position()};


        auto qdraw = [this] (const Book &book, const Reader &) -> double {
            return book.quality() + writer_quality_sd * stdnormal(Random::rng());
        };
        newbook = simulation()->create<Book>(bookPos, sharedSelf(), wrote_.size(), create_price_, create_quality_, qdraw);

        /// If enabled, add some noise in a random direction to the position
        if (writer_book_sd > 0) {
            std::normal_distribution<double> step_dist{0, writer_book_sd};
            newbook->moveBy(step_dist(Random::rng()) * Position::random(bookPos.dimensions));
        }

        wrote_.push_back(newbook);
        wrote_on_market_.insert(newbook);
    }
    ERIS_DBG("");

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

void Reader::updateBeliefs() {
    updateQualityBelief();
    updateDemandBelief();
    updateProfitStreamBelief();
    updateProfitBelief();
}

void Reader::updateQualityBelief() {
    // If we obtained any new books, update the demand belief with them
    if (not newBooks().empty()) {
        ERIS_DBG("Updating quality beliefs; prior beta = " << quality_belief_.DEBUG_betaPrior());
        std::vector<SharedMember<Book>> books;
        Eigen::VectorXd y(newBooks().size());
        size_t i = 0;
        for (auto &a : newBooks()) {
            books.push_back(a);
            y[i++] = library_[a];
        }

        quality_belief_ = quality_belief_.update(y, quality_belief_.bookData(books));
        ERIS_DBG("post beta = " << quality_belief_.DEBUG_betaPrior());
    }
}
void Reader::updateDemandBelief() {
    if (not newBooks().empty()) {
        ERIS_DBG("Updating demand beliefs; prior beta = " << demand_belief_.DEBUG_betaPrior());
        std::vector<SharedMember<Book>> books;
        VectorXd y(newBooks().size());
        MatrixXd X(newBooks().size(), demand_belief_.K());
        auto t = simulation()->t();
        size_t i = 0;
        for (auto &b : newBooks()) {
            y[i] = b->sales(t);
            X.row(i) = demand_belief_.bookRow(b, library_[b]);
        }

        demand_belief_ = demand_belief_.update(y, X);
    }
}
void Reader::updateProfitStreamBelief() {
    std::cerr << "FIXME: updateProfitStreamBelief not yet implemented!\n";
}
void Reader::updateProfitBelief() {
    std::cerr << "FIXME: updateProfitBelief not yet implemented!\n";
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
    // net utility (i.e. u-p) is highest.
    std::map<double, std::vector<SharedMember<Book>>> book_net_u;
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
        ERIS_DBG("Gonna buy a book!\n");
        // Perform some safety checks:
        // - if we've already added the book to the set of books to buy, don't add it again. (This
        //   should only be possible if there are multiple BookMarkets selling the same book).
        // - if the BookMarket says the book is no longer feasible, ignore it. (This shouldn't
        //   happen, but just in case).
        // - if we can't afford it, skip it.
        auto pinfo = bm->price(1);
        if (buy.count(book) == 0 and pinfo.feasible and pinfo.total <= money) {
            ERIS_DBG("Buying that book!");
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
    }
    reservations_.clear();
    new_books_.clear();
    std::swap(new_books_, reserved_books_); // Clear away into new_books_

    // "Eat" any money leftover
    double money = assets().remove(MONEY);
    // Finalize utility
    u_curr_ = u(money, new_books_);
    u_lifetime_ += u_curr_;
}
void Reader::intraReset() {
    auto lock = writeLock();
    reservations_.clear();
    reserved_books_.clear();
}

}
