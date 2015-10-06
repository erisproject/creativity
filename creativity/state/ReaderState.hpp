#pragma once
#include <set>
#include <map>
#include <eris/types.hpp>
#include <eris/Position.hpp>
#include "creativity/BookCopy.hpp" // IWYU pragma: keep
#include "creativity/belief/Demand.hpp"
#include "creativity/belief/Profit.hpp"
#include "creativity/belief/ProfitStream.hpp"

namespace creativity { class Reader; }

namespace creativity { namespace state {

/** Records the various variables associated with a reader.  This is basically a container class
 * with a constructor that copies the current state of a given Reader.
 */
class ReaderState final {
    public:
        /// Constructs a new ReaderState, settings its values using the given Reader.
        explicit ReaderState(const Reader &r);

        /** Constructs a new blank ReaderState for a reader with a position of the given number of
         * dimensions.  All values will be default initialized.  (The number of dimensions is needed
         * for Position initialization).
         */
        explicit ReaderState(const unsigned int dimensions);

        /// Unique simulation ID of the reader
        eris::eris_id_t id;

        /// Position of the reader
        eris::Position position;

        /** The reader's library: the keys are the book IDs of owned books, the values are the
         * per-reader specific BookCopy values.
         */
        std::map<eris::eris_id_t, BookCopy> library;

        /** The number of market-purchased books in the reader's library.  When loading/modifying a
         * ReaderState object you should call updateLibraryCounts() to recalculate this, or else set
         * it yourself. */
        unsigned int library_purchased;

        /** The number of market-purchased books in the reader's library that were acquired in the
         * current period.  When loading/modifying a ReaderState object you should call
         * updateLibraryCounts() to recalculate this, or else set it yourself. */
        unsigned int library_purchased_new;

        /** The number of public provider-purchased books in the reader's library.  When
         * loading/modifying a ReaderState object you should call updateLibraryCounts() to
         * recalculate this, or else set it yourself.
         */
        unsigned int library_public;

        /** The number of public provider-purchased books in the reader's library that were acquired
         * in the current period.  When loading/modifying a ReaderState object you should call
         * updateLibraryCounts() to recalculate this, or else set it yourself.
         */
        unsigned int library_public_new;

        /** The number of pirated books in the reader's library.  When loading/modifying a
         * ReaderState object you should call updateLibraryCounts() to recalculate this, or else set
         * it yourself. */
        unsigned int library_pirated;

        /** The number of pirated books in the reader's library that were acquired in the current
         * period.  When loading/modifying a ReaderState object you should call
         * updateLibraryCounts() to recalculate this, or else set it yourself. */
        unsigned int library_pirated_new;

        std::set<eris::eris_id_t>
            friends, ///< Friends of the reader
            new_books, ///< The set of book IDs that were newly obtained in the period, not including self-authored books.
            wrote; ///< IDs of books written by this reader

        /// Utility in the current period.
        double u;

        /// Lifetime cumulative utility up to and including the current period.
        double u_lifetime;

        /// Creation shape coefficient
        double creation_shape;

        /// Creation scale coefficient
        double creation_scale;

        // Beliefs, copied out at the time the ReaderState object is created
        belief::Profit profit; ///< Profit beliefs
        belief::Profit profit_extrap; ///< Profit beliefs using extrapolation for on-market books (will be default-constructed if there is no extrapolation)
        belief::Demand demand; ///< Single-period demand belief
        std::map<unsigned int, belief::ProfitStream> profit_stream; ///< Profit stream beliefs

        /** Returns profit_extrap if it exists (i.e. if there is extrapolated data), otherwise
         * returns profit (i.e. if there was no extrapolated data and so profit_extrap_ is exactly
         * the same as profit_).
         */
        const belief::Profit& profitExtrap() const;

        /** Utility method to recalculate `purchased`, `purchased_new`, `pirated`, and `pirated_new`
         * from the current `library`.
         *
         * \param t the time period of the state.  Needed to identify newly-acquired books (books
         * with `copy.acquired == t`
         */
        void updateLibraryCounts(eris::eris_time_t t);
};

}}
