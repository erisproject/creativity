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

Reader::Reader(const Position &pos, const Position &b1, const Position &b2)
    : WrappedPositional<agent::AssetAgent>(pos, b1, b2)
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

const std::vector<eris_id_t>& Reader::wrote() const {
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


void Reader::interApply() {
    auto &rng = Random::rng();
    // Flip a (weighted) coin to determine creativity:
    if (std::bernoulli_distribution{writer_prob}(rng)) {

        // The book is centered at the reader's position:
        Position bookPos{position()};

        // Initial price (FIXME):
        double FIXME_initial_price = 1;

        if (writer_book_sd > 0) {
            // Now add some noise to each dimension
            std::normal_distribution<double> book_noise(0, writer_book_sd);
            for (auto &x : bookPos) {
                x += book_noise(rng);
            }
        }

        // Finally we need a quality distribution.  Use chi²(1) for now.
        std::chi_squared_distribution<double> chisq1(1);
        auto quality = std::bind(chisq1, rng);

        wrote_.push_back(simulation()->create<Book>(bookPos, sharedSelf(), FIXME_initial_price, quality));
    }
}

void Reader::interAdvance() {
    // This is really crude: just move N(0,1) in a random direction 50% of the time.

    auto &rng = Random::rng();
    std::normal_distribution<double> stdnormal;
    std::normal_distribution<double> step_dist{0, 0.25};

    // Choose a random distance and a random direction
    if (std::bernoulli_distribution{0.5}(rng)) {
        const double distance = step_dist(rng);
        // Uniform hypersphere surface picking: draw a N(0,1) for each dimension, then normalize to
        // a unit distance.  See http://mathworld.wolfram.com/HyperspherePointPicking.html
        Position walk = Position::zero(position());
        double tss = 0.0;
        while (tss == 0.0) { // Extremely likely that this loop runs only once
            for (size_t i = 0; i < walk.dimensions; i++) {
                const double x = stdnormal(rng);
                walk[i] = x;
                tss += x*x;
            }
        }
        walk *= distance / sqrt(tss);
        moveBy(walk);
    }
}

double Reader::uBook(const SharedMember<Book> &b) const {
    double u = evalPolynomial(distance(b), u_poly_);
    double q_DEBUG = quality(b);
    //std::cerr << "Evaluating book " << b << ": quality=" << q_DEBUG << ", totalu=" << q_DEBUG + u << "\n";
    u += q_DEBUG;//quality(b);
    if (u < 0) u = 0.0;
    return u;
}

double Reader::quality(const SharedMember<Book> &b) const {
    auto found = library_.find(b);
    if (found != library_.end())
        return found->second;

    // Predict quality
    double e_qual = 0;
    auto author = b->author();
    auto n = author->wrote().size();
    auto coef = q_coef_.find(author);
    if (coef == q_coef_.end())
        coef = q_coef_.find(0);

    if (coef != q_coef_.end()) {
        e_qual +=
            coef->second[0] +
            coef->second[1] * n +
            coef->second[2] * n * n;
    }

    return e_qual;
}

double Reader::penalty(unsigned long n) const {
    return evalPolynomial(n, pen_poly_);
}

const std::unordered_map<eris_id_t, double>& Reader::library() { return library_; }

void Reader::intraInitialize() {
    for (auto &bm : NEW_BOOKS) {
        book_cache_.insert(bm);
    }
}

void Reader::intraOptimize() {
    auto lock = readLock();

    std::vector<SharedMember<BookMarket>> cache_del;
    // Map u-p values to sets of market ids
    std::map<double, std::vector<SharedMember<BookMarket>>> book_net_u;
    for (auto &bm : book_cache_) {
        if (not bm) { // Market removed: remove from cache and continue
            cache_del.push_back(bm);
            continue;
        }
        auto book = bm->book();
        if (library_.count(book) > 0) {
            // Already have the book: remove from cache and continue
            cache_del.push_back(bm);
            continue;
        }
        double u = uBook(book);
        double p = bm->price();
        if (u >= p) {
            // Good: the book is available and its utility (before any penalty) exceeds the purchase price
            // Store the net utility of this book:
            book_net_u[u - bm->price()].push_back(std::move(bm));
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
        auto bm = std::move(s.back());
        s.pop_back();
        // If we just pulled off the last book at the current net utility level, delete the level.
        if (s.empty()) book_net_u.erase(book_net_u.begin());

        lock.add(bm);
        auto book = bm->book();
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
        library_.emplace(book, book->qualityDraw(*this));
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
