#include "creativity/Book.hpp"
#include "creativity/BookMarket.hpp"
#include "creativity/Reader.hpp"
#include "creativity/Creativity.hpp"
#include <algorithm>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>
#include <eris/Random.hpp>


using namespace eris;

namespace creativity {

Book::Book(
        std::shared_ptr<Creativity> creativity,
        const Position &p,
        SharedMember<Reader> author,
        unsigned int order,
        double quality)
    : WrappedPositional<Good::Discrete>(p, author->wrapLowerBound(), author->wrapUpperBound()),
        creativity_{std::move(creativity)},
        author_{std::move(author)},
        order_{order},
        quality_draw_(quality, creativity_->parameters.book_quality_sd)
{}

void Book::added() {
    auto sim = simulation();
    created_ = sim->t();
    left_private_market_ = created_;

    // Make sure everything is initialized appropriately
    copies_public_total_ = copies_private_total_ = copies_pirated_total_ = 0;
    copies_sold_.clear();
    copies_pirated_.clear();
    revenue_public_total_ = revenue_private_total_ = 0;
    revenue_.clear();

    author_->registerAuthoredBook(sharedSelf());

    creativity_->newBooks().first.push_back(sharedSelf());
}

void Book::setMarket(SharedMember<BookMarket> market) {
    if (market_id_) throw std::runtime_error("Attempt to set a market for a book that already has a market");
    if (market->isPublic()) {
        market_private_ = false;
        if (public_market_created_ == 0)
            public_market_created_ = simulation()->t();
    }
    else {
        market_private_ = true;
    }
    market_id_ = market;
    dependsWeaklyOn(market);

    author_->registerMarketUpdate(sharedSelf());
}

void Book::weakDepRemoved(SharedMember<Member> mkt, eris_id_t old) {
    if (old == market_id_) {
        market_id_ = 0;
        if (not SharedMember<BookMarket>(mkt)->isPublic())
            left_private_market_ = simulation()->t();
    }

    author_->registerMarketUpdate(sharedSelf());
}

eris_time_t Book::age() const {
    return simulation()->t() - created_;
}

const eris_time_t& Book::created() const {
    return created_;
}

const eris_time_t& Book::leftPrivateMarket() const {
    return left_private_market_;
}

const eris_time_t& Book::publicMarketCreated() const {
    return public_market_created_;
}

unsigned int Book::privateMarketPeriods() const {
    auto lock = readLock();
    return hasPrivateMarket()
        ? age() + 1 // Plus 1 to count the current period
        : leftPrivateMarket() - created();
}

const unsigned int& Book::order() const {
    return order_;
}

unsigned int Book::lifeSales() const {
    auto lock = readLock();
    return lifeSalesPrivate() + lifeSalesPublic();
}

unsigned int Book::lifeSalesPrivate() const {
    return copies_private_total_;
}
unsigned int Book::lifeSalesPublic() const {
    return copies_public_total_;
}

unsigned int Book::currSales() const {
    return sales(simulation()->t());
}

unsigned int Book::sales(eris_time_t t) const {
    if (t < created_) return 0;
    auto lock = readLock();
    auto it = copies_sold_.find(t);
    if (it == copies_sold_.end()) return 0;
    return it->second;
}

eris_time_t Book::lastSale() const {
    for (auto it = copies_sold_.rbegin(); it != copies_sold_.rend(); it++) {
        if (it->second > 0) return it->first;
    }
    // No sales ever
    return 0;
}

unsigned int Book::lifePirated() const {
    auto lock = readLock();
    return copies_pirated_total_;
}
unsigned int Book::currPirated() const {
    return pirated(simulation()->t());
}
unsigned int Book::pirated(eris_time_t t) const {
    if (t < created_) return 0;
    auto lock = readLock();
    auto it = copies_pirated_.find(t);
    if (it == copies_pirated_.end()) return 0;
    return it->second;
}

void Book::recordSale(unsigned int count, double revenue) {
    auto lock = writeLock();
    if (not hasAnyMarket())
        throw std::runtime_error("Book: cannot record sales for a book without a market");

    if (hasPrivateMarket()) {
        copies_private_total_ += count;
        revenue_private_total_ += revenue;
    }
    else {
        copies_public_total_ += count;
        revenue_public_total_ += revenue;
    }
    auto t = simulation()->t();
    copies_sold_[t] += count;
    revenue_[t] += revenue;
}

void Book::recordPiracy(unsigned int new_copies) {
    auto lock = writeLock();
    copies_pirated_total_ += new_copies;
    copies_pirated_[simulation()->t()] += new_copies;
}

void Book::recordPrize(double prize) {
    auto lock = writeLock();

    prize_[simulation()->t()] += prize;
    prize_total_ += prize;
}

double Book::lifeRevenuePrivate() const {
    return revenue_private_total_;
}

double Book::lifeRevenuePublic() const {
    return revenue_public_total_;
}

double Book::lifeRevenue() const {
    auto lock = readLock();
    return revenue_public_total_ + revenue_private_total_;
}

double Book::lifeProfitPrivate() const {
    return lifeRevenuePrivate()
        - privateMarketPeriods() * creativity_->parameters.cost_market
        - lifeSalesPrivate() * creativity_->parameters.cost_unit;
}

double Book::currRevenue() const {
    return revenue(simulation()->t());
}

double Book::revenue(eris_time_t t) const {
    if (t < created_) return 0.0;
    auto lock = readLock();
    auto it = revenue_.find(t);
    if (it == revenue_.end()) return 0.0;
    return it->second;
}

double Book::currPrize() const {
    return prize(simulation()->t());
}

double Book::prize(eris_time_t t) const {
    if (t < created_) return 0.0;
    auto lock = readLock();
    auto it = prize_.find(t);
    if (it == prize_.end()) return 0.0;
    return it->second;
}

double Book::lifePrize() const {
    return prize_total_;
}

unsigned int Book::copies() const {
    auto lock = readLock();
    return lifeSales() + lifePirated();
}

unsigned int Book::queryCopies() const {
    unsigned int copies = 0;
    SharedMember<Book> me(sharedSelf());
    for (auto &r : simulation()->agents<Reader>()) {
        if (r->library().count(me))
            copies++;
    }
    return copies;
}

bool Book::livingAuthor() const {
    return author_->hasSimulation();
}

SharedMember<Reader> Book::author() const {
    return author_;
}

bool Book::hasAnyMarket() const {
    return market_id_ != 0;
}

bool Book::hasPrivateMarket() const {
    auto lock = readLock();
    return market_private_ and market_id_ != 0;
}

bool Book::hasPublicMarket() const {
    auto lock = readLock();
    return not market_private_ and market_id_ != 0;
}

SharedMember<BookMarket> Book::market() const {
    return simMarket<BookMarket>(market_id_);
}

double Book::price() const {
    return hasAnyMarket()
        ? market()->price()
        : std::numeric_limits<double>::quiet_NaN();
}

double Book::qualityMean() const {
    return quality_draw_.mean();
}

double Book::qualityDraw() {
    // We could call Random::truncDist here, but with its defaults for a normal it's just going to
    // end up doing a rejection sampling loop anyway, so just do that directly.  This loop shouldn't
    // run too long; since the lowest possible quality is 0, the worst case is that this draw
    // produces invalid values half the time.
    if (quality_draw_.stddev() > 0) {
        auto &rng = Random::rng();
        double x;
        do { x = quality_draw_(rng); } while (x < 0);
        return x;
    }
    // No stddev: just return the exact quality value
    return quality_draw_.mean();
}

}
