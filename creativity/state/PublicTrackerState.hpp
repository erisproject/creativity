#pragma once
#include <eris/types.hpp>

namespace creativity { class PublicTracker; }

namespace creativity { namespace state {

/** Records the various variables associated with a reader.  This is basically a container class
 * with a constructor that copies the current state of a given Reader.
 */
class PublicTrackerState final {
    public:
        /// Creates a state representing the given PublicTracker
        PublicTrackerState(const PublicTracker &pt);

        /// Default constructor: creates a PublicTrackerState with default-initialized fields.
        PublicTrackerState() = default;

        /// Unique simulation ID of the reader
        eris::eris_id_t id;

        /// The per-period lump-sum tax amount for download-count public sharing
        double dl_tax;

        /// The per-period lump-sum tax amount for vote-count public sharing
        double vote_tax;

        /** Assets unspent at the end of the period.  This will be non-zero at the end of a period
         * only if the period had no public tracker sales, and thus no authors to distribute the tax
         * to (in which case the leftover assets are added to the next period's collected assets).
         */
        double dl_unspent;
        
        /** Assets unspent at the end of the period for the vote-based redistribution.  Like
         * `dl_unspent`, this should only be non-zero when enabled and there are no votes in a
         * period.
         */
        double vote_unspent;
};

}}
