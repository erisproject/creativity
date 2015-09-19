#include "creativity/BookMarket.hpp"
#include "creativity/Book.hpp"
#include "creativity/Reader.hpp"
#include "creativity/Creativity.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>

using eris::agent::AssetAgent;
using namespace eris;

namespace creativity {

BookMarket::BookMarket(std::shared_ptr<Creativity> creativity, SharedMember<Book> b, double price)
    : Market{{{ b, 1 }}, {{ creativity->money, 1 }}}, creativity_{std::move(creativity)}, book_{b}, price_{price}
{}

SharedMember<Book> BookMarket::book() {
    return book_;
}

void BookMarket::added() {
    book_->setMarket(sharedSelf());
}

BookMarket::price_info BookMarket::price(double q) const {
    unsigned long ql = (unsigned long) std::floor(q);
    return price_info(ql*price_, price_, price_);
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

void BookMarket::buy(Reservation &res) {
    if (res.state != ReservationState::pending)
        throw Reservation::non_pending_exception();

    auto lock = writeLock(book_);
    Bundle &b = reservationBundle_(res);

    // Record the sale in the book status
    book_->recordSale(res.quantity, b[creativity_->money]);

    // Transfer the money into the "proceeds" jar (which will eventually go to the author)
    b.transferApprox(b, proceeds_);

    // Leave a copy of the book
    b.set(book_, res.quantity);

    // Complete the reservation
    Market::buy(res);
}

void BookMarket::intraFinish() {
    auto author = book_->author();
    auto lock = writeLock(book_, author);
    // First subtract off variable costs incurred
    BundleNegative tvc(creativity_->money, book_->currSales() * -creativity_->parameters.cost_unit);
    tvc.transferApprox(tvc, proceeds_, 1e-8);

    // Transfer profits to the author
    author->assets() += proceeds_;

    proceeds_.clear();
}

}
