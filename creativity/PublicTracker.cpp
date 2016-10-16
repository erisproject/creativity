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
    bool dl = creativity_.parameters.policy & POLICY_PUBLIC_SHARING;
    bool vote = creativity_.parameters.policy & POLICY_PUBLIC_SHARING_VOTING;
    if (!dl && !vote)
        throw std::logic_error("PublicTracker created for a simulation without public sharing or public voting!");

    if (dl and dlTax() < 0) throw std::domain_error("PublicTracker creation error: public sharing lump sum tax cannot be negative");
    if (vote) {
        if (voteTax() <= 0) throw std::domain_error("PublicTracker creation error: public voting lump sum tax must be positive");
        if (creativity_.parameters.policy_public_sharing_voting_votes <= 0)
            throw std::domain_error("PublicTracker creation error: number of votes must be positive");
    }
}

void PublicTracker::interApply() {
    auto lock = writeLock();
    // If the tax is 0, we don't need to do any transfers.
    if (creativity_.publicSharing() and dlTax() > 0) {
        Bundle tax_bill(creativity_.money, dlTax());
        for (auto &r : simulation()->agents<Reader>()) {
            lock.add(r);
            r->assets.transferApprox(tax_bill, dl_assets, 1e-6);
            lock.remove(r);
        }
    }

    if (creativity_.publicSharingVoting() and voteTax() > 0) {
        Bundle tax_bill(creativity_.money, voteTax());
        for (auto &r : simulation()->agents<Reader>()) {
            lock.add(r);
            r->assets.transferApprox(tax_bill, vote_assets, 1e-6);
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

void PublicTracker::distributeDLFunds() {
    unsigned int total_copies = 0;
    std::unordered_map<SharedMember<Reader>, uint32_t> author_copies;
    std::unordered_map<SharedMember<Book>, uint32_t> book_copies;
    if (dl_assets[creativity_.money] > 0) {
        // Get the number of PublicTrackerMarket sales for each author
        for (auto &ptm : simulation()->markets<PublicTrackerMarket>()) {
            auto copies = ptm->book()->currSales();
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

void PublicTracker::distributeVoteFunds() {
    unsigned int total_votes = 0;
    std::unordered_map<SharedMember<Reader>, uint32_t> author_votes;
    std::unordered_map<SharedMember<Book>, uint32_t> book_votes;
    if (vote_assets[creativity_.money] > 0) {
        // Get the number of PublicTrackerMarket votes for each author
        for (auto &ptm : simulation()->markets<PublicTrackerMarket>()) {
            auto votes = ptm->book()->currVotes();
            if (votes > 0) {
                total_votes += votes;
                author_votes[ptm->book()->author()] += votes;
                book_votes[ptm->book()] += votes;
            }
        }
    }

    Bundle per_vote_payout(creativity_.money, assets[creativity_.money] / total_votes);

    if (total_votes > 0) {
        for (auto &ac : author_votes) {
            assets.transferApprox(per_vote_payout * ac.second, ac.first->assets, 1e-6);
        }
        for (auto &bc : book_votes) {
            bc.first->recordPrize(per_vote_payout[creativity_.money] * bc.second);
        }
    }
    // Otherwise no votes at all: just leave the assets for the next period's pool
}

void PublicTracker::intraFinish() {

    if (creativity_.publicSharing())
        distributeDLFunds();

    if (creativity_.publicSharingVoting())
        distributeVoteFunds();

}

const double& PublicTracker::votePool() const {
    return vote_assets[creativity_.money];
}

const double& PublicTracker::dlPool() const {
    return dl_assets[creativity_.money];
}

}
