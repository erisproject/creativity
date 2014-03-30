#include "creativity/Book.hpp"
#include "creativity/BookMarket.hpp"
#include "creativity/Reader.hpp"

using namespace eris;

namespace creativity {

Book::Book(
        Position p,
        SharedMember<Reader> author,
        double initial_price,
        std::function<double(const Book&, const Reader&)> quality)
    : Positional<Good::Discrete>(p),
        author_{author},
        init_price_{initial_price},
        quality_{quality}
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
    if (old == author_) author_ = 0;
    else if (old == market_) market_ = 0;
}

unsigned long Book::age() const {
    return simulation()->t() - created_;
}

unsigned long Book::sales() const {
    auto lock = readLock();
    return copies_sold_;
}

void Book::sales(unsigned long new_sales) {
    auto lock = writeLock();
    copies_sold_ += new_sales;
}

bool Book::hasAuthor() const {
    return author_ != 0;
}

eris_id_t Book::authorID() const {
    return author_;
}

SharedMember<Reader> Book::author() const {
    return simAgent<Reader>(author_);
}

bool Book::onMarket() const {
    return market_ != 0;
}

SharedMember<BookMarket> Book::market() const {
    return simMarket<BookMarket>(market_);
}

double Book::qualityDraw(const Reader &reader) {
    return std::max(0.0, quality_(*this, reader));
}

}
