#include "creativity/Book.hpp"
#include "creativity/BookMarket.hpp"
#include "creativity/Reader.hpp"
#include "creativity/Creativity.hpp"

using namespace eris;

namespace creativity {

Book::Book(
        std::shared_ptr<Creativity> creativity,
        const Position &p,
        SharedMember<Reader> author,
        unsigned int order,
        double quality,
        std::function<double(const Book&)> qDraw)
    : WrappedPositional<Good::Discrete>(p, author->wrapLowerBound(), author->wrapUpperBound()),
        creativity_{std::move(creativity)},
        author_{std::move(author)},
        order_{order},
        quality_{quality},
        quality_draw_{std::move(qDraw)}
{}

void Book::added() {
    auto sim = simulation();
    created_ = sim->t();

    // These shouldn't be doing anything, but clear everything just in case a Book is removed and reintroduced:
    copies_sold_total_ = 0;
    copies_sold_.clear();
    copies_pirated_total_ = 0;
    copies_pirated_.clear();
    revenue_total_ = 0;
    revenue_.clear();

    creativity_->newBooks().first.push_back(sharedSelf());
}

void Book::setMarket(SharedMember<BookMarket> market, bool primary) {
    if (market_) throw std::runtime_error("Attempt to set a market for a book that already has a market");
    market_primary_ = primary;
    market_ = market;
    out_of_print_ = 0;
    dependsWeaklyOn(market);
}

void Book::weakDepRemoved(SharedMember<Member>, eris_id_t old) {
    if (old == market_) {
        market_ = 0;
        out_of_print_ = simulation()->t();
    }
}

eris_time_t Book::age() const {
    return simulation()->t() - created_;
}

const eris_time_t& Book::created() const {
    return created_;
}

const eris_time_t& Book::outOfPrint() const {
    return out_of_print_;
}

unsigned int Book::marketPeriods() const {
    return hasMarket()
        ? age() + 1 // Plus 1 to count the current period
        : outOfPrint() - created();
}

const unsigned int& Book::order() const {
    return order_;
}

unsigned int Book::lifeSales() const {
    auto lock = readLock();
    return copies_sold_total_;
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
    copies_sold_total_ += count;
    revenue_total_ += revenue;
    auto t = simulation()->t();
    copies_sold_[t] += count;
    revenue_[t] += revenue;
}

void Book::recordPiracy(unsigned int new_copies) {
    auto lock = writeLock();
    copies_pirated_total_ += new_copies;
    copies_pirated_[simulation()->t()] += new_copies;
}

double Book::lifeRevenue() const {
    auto lock = readLock();
    return revenue_total_;
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

bool Book::hasMarket() const {
    return market_ != 0;
}

bool Book::hasPrimaryMarket() const {
    return market_primary_ and market_ != 0;
}

SharedMember<BookMarket> Book::market() const {
    return simMarket<BookMarket>(market_);
}

double Book::price() const {
    return hasMarket()
        ? market()->price()
        : std::numeric_limits<double>::quiet_NaN();
}

const double& Book::quality() const {
    return quality_;
}

double Book::qualityDraw() {
    return std::max(0.0, quality_draw_(*this));
}

}
