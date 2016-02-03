#include "creativity/PublicTrackerMarket.hpp"
#include "creativity/Creativity.hpp"
#include "creativity/Book.hpp"
#include <algorithm>

namespace creativity {

using namespace eris;

PublicTrackerMarket::PublicTrackerMarket(const Creativity &creativity, SharedMember<Book> b)
    : BookMarket(creativity, b, 0.0) {}

void PublicTrackerMarket::added() {
    BookMarket::added();
    updatePrice();
}

void PublicTrackerMarket::updatePrice() {
    setPrice(std::min(creativity_.parameters.cost_unit, creativity_.parameters.cost_piracy));
}

void PublicTrackerMarket::interAdvance() {
    updatePrice();
}

void PublicTrackerMarket::intraFinish() {
    // This is overriding BookMarket to skip everything it does: price equals marginal cost, and so
    // there can never be any profits to transfer.  Clear the proceeds_ jar (in case it contains
    // epsilon amounts).
    proceeds_.clear();
}

}
