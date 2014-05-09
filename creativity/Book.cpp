#include "creativity/Book.hpp"
#include "creativity/BookMarket.hpp"
#include "creativity/Reader.hpp"

using namespace eris;

namespace creativity {

Book::Book(
        Position p,
        SharedMember<Reader> author,
        unsigned long order,
        double initial_price,
        double quality,
        std::function<double(const Book&, const Reader&)> qDraw)
    : WrappedPositional<Good::Discrete>(p, author->wrapLowerBound(), author->wrapUpperBound()),
        author_{std::move(author)},
        order_{std::move(order)},
        init_price_{std::move(initial_price)},
        quality_{std::move(quality)},
        quality_draw_{std::move(qDraw)}
{}

void Book::added() {
    auto sim = simulation();
    created_ = sim->t();
    copies_sold_ = 0;
    auto mkt = sim->create<BookMarket>(*this, init_price_);
    NEW_BOOKS.push_back(mkt);
    market_ = mkt;
    dependsWeaklyOn(mkt);
}

void Book::weakDepRemoved(SharedMember<Member>, const eris::eris_id_t &old) {
    if (old == market_) market_ = 0;
}

unsigned long Book::age() const {
    return simulation()->t() - created_;
}

const unsigned long& Book::order() const {
    return order_;
}

unsigned long Book::sales() const {
    auto lock = readLock();
    return copies_sold_;
}

void Book::sales(unsigned long new_sales) {
    auto lock = writeLock();
    copies_sold_ += new_sales;
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

const double& Book::quality() const {
    return quality_;
}

double Book::qualityDraw(const Reader &reader) {
    return std::max(0.0, quality_draw_(*this, reader));
}

}
