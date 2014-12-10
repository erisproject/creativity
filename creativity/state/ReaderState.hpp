#pragma once
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <eris/types.hpp>
#include <eris/Position.hpp>
#include <eris/noncopyable.hpp>
#include "creativity/Reader.hpp"

namespace creativity { namespace state {

/** Records the various variables associated with a reader.  This is basically a container class
 * with a constructor that copies the current state of a given Reader.
 */
class ReaderState final {
    public:
        /// Constructs a new ReaderState, settings its values using the given Reader.
        ReaderState(const Reader &r);

        /** Constructs a new blank ReaderState for a reader with a position of the given number of
         * dimensions.  All values will be default initialized.  (The number of dimensions is needed
         * for Position initialization).
         */
        ReaderState(const unsigned int dimensions);

        /// Unique simulation ID of the reader
        eris::eris_id_t id;

        /// Position of the reader
        eris::Position position;

        /** The reader's library: the keys are the book IDs of owned books, the values are the book
         * quality values realized by this reader.
         */
        std::unordered_map<eris::eris_id_t, double> library;

        std::unordered_set<eris::eris_id_t>
            friends, ///< Friends of the reader
            library_purchased, ///< The set of book IDs in `library` that were purchased.
            library_pirated, ///< The set of book IDs in `library` that were pirated.
            new_books, ///< The set of book IDs that were newly obtained in the period, not including self-authored books.
            new_purchased, ///< The set of new book IDs that were purchased in the period.
            new_pirated; ///< The set of new book IDs that were pirated in the period.

        /** Set of IDs of books written by this reader, sorted by ID (and thus also by creation
         * order). */
        std::set<eris::eris_id_t> wrote;

        /// Utility in the current period.
        double u;

        /// Lifetime cumulative utility up to and including the current period.
        double u_lifetime;

        /// Fixed cost of keeping a book on the market
        double cost_fixed;

        /// Unit cost of creating a copy of a book
        double cost_unit;

        /// Unit cost of obtaining a pirated copy of a book
        double cost_piracy;

        /// Per-period (external) income
        double income;

        // Beliefs, copied out at the time the ReaderState object is created
        belief::Profit profit; ///< Profit beliefs
        belief::Profit profit_extrap; ///< Profit beliefs using extrapolation for on-market books
        belief::Demand demand; ///< Single-period demand belief
        belief::Quality quality; ///< Quality belief
        std::map<unsigned int, belief::ProfitStream> profit_stream; ///< Profit stream beliefs

    private:
        // Copy from one container of SharedMembers into a container of IDs
        template <class F, class T>
        void copyIDs(const F &from, T &to) {
            to.reserve(from.size());
            for (const auto &m : from) to.emplace(m->id());
        }
};

}}
