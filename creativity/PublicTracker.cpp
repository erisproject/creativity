#include "creativity/PublicTracker.hpp"
#include "creativity/PublicTrackerMarket.hpp"
#include "creativity/Reader.hpp"
#include "creativity/Book.hpp"
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace creativity {

using namespace eris;

PublicTracker::PublicTracker(const Creativity &creativity) : creativity_{std::move(creativity)} {
    if (tax() < 0) throw std::domain_error("PublicTracker creation error: lump sum tax cannot be negative");
}

void PublicTracker::interApply() {
    auto lock = writeLock();
    // If the tax is 0, we don't need to do any transfers.
    if (tax() > 0) {
        Bundle tax_bill(creativity_.money, tax());
        for (auto &r : simulation()->agents<Reader>()) {
            lock.add(r);
            r->assets.transferApprox(tax_bill, assets, 1e-6);
            lock.remove(r);
        }
    }

    // Create new markets for off-market books
    for (auto &b : simulation()->goods<Book>()) {
        if (not b->hasAnyMarket()) { // The author decided not to put the book on the market for the upcoming period, so we'll take over
            simulation()->spawn<PublicTrackerMarket>(creativity_, b);
        }
    }
}

void PublicTracker::intraFinish() {
    unsigned int total_copies = 0;
    std::unordered_map<SharedMember<Reader>, unsigned int> author_copies;
    std::unordered_map<SharedMember<Book>, unsigned int> book_copies;
    if (assets[creativity_.money] > 0) {
        // Get the number of PublicTrackerMarket sales for each author
        for (auto &ptm : simulation()->markets<PublicTrackerMarket>()) {
            unsigned int copies = ptm->book()->currSales();
            if (copies > 0) {
                total_copies += copies;
                author_copies[ptm->book()->author()] += copies;
                book_copies[ptm->book()] += copies;
            }
        }
    }

    Bundle per_copy_payout(creativity_.money, assets[creativity_.money] / total_copies);

    if (total_copies > 0) {
        for (auto &ac : author_copies) {
            assets.transferApprox(per_copy_payout * ac.second, ac.first->assets, 1e-6);
        }
        for (auto &bc : book_copies) {
            bc.first->recordPrize(per_copy_payout[creativity_.money] * bc.second);
        }
    }
    // Otherwise no copies downloaded: just leave the assets for the next period's pool
}

const double& PublicTracker::pool() const {
    return assets[creativity_.money];
}

}
