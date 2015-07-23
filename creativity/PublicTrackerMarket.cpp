#include "creativity/PublicTrackerMarket.hpp"
#include "creativity/Creativity.hpp"
#include "creativity/CreativitySettings.hpp"
#include <algorithm>

namespace creativity {

using namespace eris;

PublicTrackerMarket::PublicTrackerMarket(std::shared_ptr<Creativity> creativity, SharedMember<Book> b)
    : BookMarket(creativity, b, 0.0) {
    updatePrice();
}

void PublicTrackerMarket::updatePrice() {
    setPrice(std::min(creativity_->parameters.cost_unit, creativity_->parameters.cost_piracy));
}

void PublicTrackerMarket::interAdvance() {
    updatePrice();
}

void PublicTrackerMarket::intraFinish() {
    proceeds_.clear();
}

}
