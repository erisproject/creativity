#pragma once
#include <eris/Agent.hpp>
#include <eris/Optimize.hpp>
#include "creativity/Creativity.hpp"

namespace creativity {

class Creativity;

/** This class represents a public provider of book copies that pays author based on the number of
 * downloads or votes received in a period.  In particular, it collects a lump sum tax amount just
 * after each agent receives his non-book income, then pays out this amount at the end of the period
 * to authors whose books were downloaded, proportional to the number of downloads.
 *
 * Per-unit costs are whichever of `cost_unit` or `cost_piracy` is lowest (that is, the public
 * provider is always able to use the least-cost means of distribution), and copies are priced at
 * this marginal cost.
 *
 * The lump sum tax can be 0, in which case out-of-market books are provided at marginal cost
 * without any payment going to authors.
 *
 * For example, if there are 100 agents, the tax amount is 10 (so 1000 in total), and 4 books (A-D)
 * are downloaded in the period according to:
 *     A - 15 downloads
 *     B - 5 downloads
 *     C - 4 downloads
 *     D - 1 download
 * then A's author would receive 600, B's author 200, C's author 160, and D's author 40.
 *
 * If there are no public copies downloaded in a period at all, the lump sum amount is retained for
 * the next period: thus any authors downloaded in the next period will receive a double payoff.
 *
 * When operating in vote mode the payoff is the same, but is based on the number of votes rather
 * than the number of downloads.
 *
 * Note that if both public sharing *and* public sharing with voting are enabled, both lump sum
 * taxes are collected, pooled separately, and redistributed according to each method.
 *
 * This tracker automatically takes over any books that leave the market; authors may still choose
 * to keep books on the market themselves, but doing so risks piracy, and receives no compensation
 * from this public provider.
 *
 * Authors may also, when this provider is available, choose not to put books on the market at all
 * but to write and release directly to the public tracker.
 */
class PublicTracker : public eris::Agent,
    public virtual eris::interopt::Apply, public virtual eris::intraopt::Finish {
    public:
        PublicTracker() = delete; ///< Not default constructible

        /** Constructor for a new PublicTracker agent.  Takes a reference to the creativity object.
         *
         * \throws std::domain_error if `creativity.policy_public_sharing_tax` is negative.
         * \throws std::logic_error if `creativity.policy` has neither public sharing nor public
         * sharing with voting policies enabled.
         */
        explicit PublicTracker(const Creativity &creativity);

        /** When the period advances, we take the lump sum tax from all agents, and create a public
         * market for any books that don't have a market (i.e. the author isn't marketing them
         * anymore). */
        void interApply() override;

        /** Override priority to run after the Reader's interApply has deposited income, and after
         * readers have created new books and withdrawn old books from the market.
         */
        double interApplyPriority() const override { return 1.0; }

        /** When the period finishes, we return the lump sum tax proportionally to all authors. */
        void intraFinish() override;

        /// Returns the lump sum per-reader tax collected each period for download count compensation
        double dlTax() const { return creativity_.parameters.policy_public_sharing_tax; }

        /// Returns the lump sum per-reader tax collected each period for vote-based compensation
        double voteTax() const { return creativity_.parameters.policy_public_voting_tax; }

        /** Returns the current total asset pool that will be distributed to authors (proportional
         * to downloads) at the end of the current period.  If called between periods, this will be
         * 0 unless the previous period had no public downloads at all: in such a case, this is the
         * amount that will be added to collected lump sum amounts in the upcoming period.
         */
        const double& dlPool() const;

        /// Like pool(), but for the voting version of the public tracker.
        const double& votePool() const;

    protected:
        /// The Creativity object that owns the simulation this reader belongs to
        const Creativity &creativity_;

        /// Collected revenue for download-count redistribution
        eris::Bundle dl_assets;

        /// Collected revenue for vote-count redistribution
        eris::Bundle vote_assets;

        /** Awards all collected download funds to authors in proportion to the number of public
         * downloads of copies of their works.
         *
         * This is called by intraFinish() when the public tracker is doing per-download awards.
         */
        void distributeDLFunds();

        /** Awards all collected voting funds to authors in proportion to the number of votes
         * received for their works.
         *
         * This is called by intraFinish() when the public tracker with voting is enabled.
         */
        void distributeVoteFunds();
};

}
