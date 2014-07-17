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
 * with a constructor that copies the current state of a given Reader. */
class ReaderState final : private eris::noncopyable {
    public:
        /** Constructs a new ReaderState without setting any of its values (they will be default
         * initialized).
         */
        ReaderState() = default;

        /// Constructs a new ReaderState, settings its values using the given Reader.
        ReaderState(const Reader &r);

        /// Unique simulation ID of the reader
        eris::eris_id_t id;

        /// The simulation period this state represents.
        unsigned long t;

        /// Position of the reader
        eris::Position position;

        /** The reader's library: the keys are the book IDs of owned books, the values are the book
         * quality values realized by this reader.
         */
        std::unordered_map<eris::eris_id_t, double> library;

        /** A set of book IDs that were newly obtained in the given period. */
        std::unordered_set<eris::eris_id_t> newBooks;

        /** The list of book id's of books that were written by this author, in order from oldest to
         * newest.
         */
        std::vector<eris::eris_id_t> wrote;

        /// Utility in the current period.
        double u;

        /// Lifetime cumulative utility up to and including the current period.
        double uLifetime;

        /// Agent beliefs, copied out at the time the ReaderState object is created
        struct {
            belief::Profit profit; ///< Profit beliefs
            belief::Profit profit_extrap; ///< Profit beliefs using extrapolation for on-market books
            belief::Demand demand; ///< Single-period demand belief
            belief::Quality quality; ///< Quality belief
            std::map<unsigned long, belief::ProfitStream> profit_stream; ///< Profit stream beliefs
        } belief;
};

}}
