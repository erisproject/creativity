#include "creativity/Book.hpp"
#include "creativity/BookMarket.hpp"

using namespace eris;

namespace creativity {

Book::Book(Position p, SharedMember<Reader> author, double initial_price)
    : Positional<Good::Discrete>(p), author_(author), init_price_{initial_price}
{}

void Book::added() {
    auto sim = simulation();
    created_ = sim->t();
    copies_sold_ = 0;
    market_ = sim->create<BookMarket>(*this, init_price_);
    dependsWeaklyOn(market_);
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

SharedMember<Reader> Book::author() const {
    return simAgent<Reader>(author_);
}

bool Book::onMarket() const {
    return market_ != 0;
}

SharedMember<BookMarket> Book::market() const {
    return simMarket<BookMarket>(market_);
}

}
