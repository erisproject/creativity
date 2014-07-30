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
        unsigned long order,
        double initial_price,
        double quality,
        std::function<double(const Book&, const Reader&)> qDraw)
    : WrappedPositional<Good::Discrete>(p, author->wrapLowerBound(), author->wrapUpperBound()),
        creativity_{std::move(creativity)},
        author_{std::move(author)},
        order_{order},
        init_price_{initial_price},
        quality_{quality},
        quality_draw_{std::move(qDraw)}
{}

void Book::added() {
    auto sim = simulation();
    created_ = sim->t();

    // These shouldn't be doing anything, but clear everything just in case a Book is removed and reintroduced:
    copies_sold_total_ = 0;
    copies_sold_.clear();
    revenue_total_ = 0;
    revenue_.clear();

    auto mkt = sim->create<BookMarket>(creativity_, sharedSelf(), init_price_);
    creativity_->newBooks().first.push_back(sharedSelf());
    market_ = mkt;
    dependsWeaklyOn(mkt);
}

void Book::weakDepRemoved(SharedMember<Member>, eris_id_t old) {
    if (old == market_) {
        market_ = 0;
        out_of_print_ = simulation()->t();
    }
}

unsigned long Book::age() const {
    return simulation()->t() - created_;
}

const unsigned long& Book::created() const {
    return created_;
}

const unsigned long& Book::outOfPrint() const {
    return out_of_print_;
}

unsigned long Book::marketPeriods() const {
    return hasMarket()
        ? age() + 1 // Plus 1 to count the current period
        : outOfPrint() - created();
}

const unsigned long& Book::order() const {
    return order_;
}

unsigned long Book::lifeSales() const {
    auto lock = readLock();
    return copies_sold_total_;
}

unsigned long Book::currSales() const {
    return sales(simulation()->t());
}

unsigned long Book::sales(unsigned long t) const {
    if (t < created_) return 0;
    auto lock = readLock();
    auto it = copies_sold_.find(t);
    if (it == copies_sold_.end()) return 0;
    return it->second;
}

void Book::sale(unsigned long count, double revenue) {
    auto lock = writeLock();
    copies_sold_total_ += count;
    revenue_total_ += revenue;
    auto t = simulation()->t();
    copies_sold_[t] += count;
    revenue_[t] += revenue;
}

double Book::lifeRevenue() const {
    auto lock = readLock();
    return revenue_total_;
}

double Book::currRevenue() const {
    return revenue(simulation()->t());
}

double Book::revenue(unsigned long t) const {
    if (t < created_) return 0.0;
    auto lock = readLock();
    auto it = revenue_.find(t);
    if (it == revenue_.end()) return 0.0;
    return it->second;
}

unsigned long Book::copies() const {
    unsigned long copies = 0;
    SharedMember<Book> me{sharedSelf()};
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

double Book::qualityDraw(const Reader &reader) {
    return std::max(0.0, quality_draw_(*this, reader));
}

}
