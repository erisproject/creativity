#include "creativity/BookMarket.hpp"
#include "creativity/Book.hpp"
#include "creativity/Reader.hpp"
#include "creativity/Creativity.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>

using namespace eris;

namespace creativity {

// NB: b might not actually be in the simulation yet (in which case it won't have an eris_id)
BookMarket::BookMarket(const Creativity &creativity, SharedMember<Book> b, double price)
    // Our output bundle is empty: book copies aren't Goods, and are handled elsewhere
    : Market{Bundle{}, {{ creativity.money, 1 }}}, creativity_{creativity}, book_{b}, price_{price}
{}

SharedMember<Book> BookMarket::book() {
    return book_;
}

void BookMarket::added() {
    book_->setMarket(sharedSelf());
}

BookMarket::price_info BookMarket::price(double q) const {
    if (q != 1) throw std::logic_error("BookMarket price() can only b called with q=1");
    return price_info(price_, price_, price_);
}

void BookMarket::setPrice(double p) {
    if (p <= 0) throw std::domain_error("Cannot call BookMarket::setPrice() with a non-positive price");
    if (hasSimulation() and simulation()->runStageIntra()) throw std::logic_error("Cannot change book price during intra-optimization");
    price_ = p;
}

const double& BookMarket::price() {
    return price_;
}

BookMarket::quantity_info BookMarket::quantity(double p) const {
    return {
        .quantity = (p >= price_ ? 1. : 0.),
        .constrained = false,
        .spent = (p >= price_ ? price_ : 0.),
        .unspent = (p >= price_ ? price_ - p : p)
    };
}

Market::Reservation BookMarket::reserve(
        SharedMember<Agent> agent, double q, double p_max) {
    if (q != 1) throw std::logic_error("BookMarket price() can only b called with q=1");

    if (price_ > p_max) throw std::logic_error("BookMarket::reserve() called with p_max < price");

    auto lock = writeLock(agent); // Lock the market and agent

    return createReservation(agent, 1.0, price_);
}

void BookMarket::buy(Reservation &res) {
    if (res.state != ReservationState::pending)
        throw Reservation::non_pending_exception();

    auto lock = writeLock(book_);
    Bundle &b = reservationBundle_(res);

    // Record the sale in the book status
    book_->recordSale(res.quantity, b[creativity_.money]);

    // Transfer the money into the "proceeds" jar (which will eventually go to the author)
    b.transferApprox(b, proceeds_);

    // We don't transfer anything back: the new BookCopy is actually created in Reader

    // Complete the reservation
    Market::buy(res);
}

void BookMarket::intraFinish() {
    auto author = book_->author();
    auto lock = writeLock(book_, author);
    // First subtract off variable costs incurred
    BundleNegative tvc(creativity_.money, book_->currSales() * -creativity_.parameters.cost_unit);
    tvc.transferApprox(tvc, proceeds_, 1e-6);

    // Transfer profits to the author
    author->assets += proceeds_;

    proceeds_.clear();
}

}
