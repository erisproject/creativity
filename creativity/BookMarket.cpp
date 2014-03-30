#include "creativity/BookMarket.hpp"
#include "creativity/Reader.hpp"
#include "creativity/Book.hpp"

using namespace eris;

namespace creativity {

BookMarket::BookMarket(const Book &b, const double &price)
    : Market{{{ b, 1 }}, {{ MONEY, 1 }}}, book_{b}, price_{price}
{}

SharedMember<Book> BookMarket::book() {
    return simGood<Book>(book_);
}

BookMarket::price_info BookMarket::price(double q) const {
    unsigned long ql = (unsigned long) std::floor(q);
    return {
        .feasible = true,
            .total = ql*price_,
            .marginal = price_,
            .marginalFirst = price_
    };
}

void BookMarket::setPrice(double p) {
    if (p <= 0) throw std::domain_error("Cannot call BookMarket::setPrice() with a non-positive price");
    if (simulation()->runStageIntra()) throw std::logic_error("Cannot change book price during intra-optimization");
    price_ = p;
}

const double& BookMarket::price() {
    return price_;
}

BookMarket::quantity_info BookMarket::quantity(double p) const {
    unsigned long ql = std::floor(p / price_);
    p -= ql*price_;
    return {
        .quantity = (double) ql,
            .constrained = false,
            .spent = ql*price_,
            .unspent = p
    };
}

Market::Reservation BookMarket::reserve(
        SharedMember<AssetAgent> agent, double q, double p_max) {

    auto ql = (unsigned long) std::floor(q);
    if (ql * price_ > p_max) {
        ql = std::floor(p_max / price_);
    }
    double total_p = ql*price_;

    auto lock = writeLock(agent); // Lock the market and agent

    return createReservation(agent, ql, total_p);
}

void BookMarket::buy_(Reservation_ &res) {
    if (res.state != ReservationState::pending)
        throw Reservation_::non_pending_exception();

    auto book = simGood<Book>(book_);
    auto lock = writeLock(book);
    Bundle &b = reservationBundle_(res);
    book->sales(res.quantity);
    b.transferApprox(b, proceeds_); // Take the money
    b.set(book_, res.quantity); // Leave a copy of the book
    Market::buy_(res);
}

void BookMarket::interAdvance() {
    auto author = simGood<Book>(book_)->author();
    auto lock = writeLock(author);
    proceeds_.transferApprox(proceeds_, author->assets());
}

}
