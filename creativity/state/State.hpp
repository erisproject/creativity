#pragma once
#include <eris/Simulation.hpp>
#include "creativity/state/ReaderState.hpp"
#include "creativity/state/BookState.hpp"
#include "creativity/state/PublicTrackerState.hpp"
#include <cstdint>
#include <limits>
#include <map>

namespace creativity { namespace state {

/** Class storing the state of the simulation at the end of a simulation period. */
class State final {
    public:
        /** Copies the current simulation readers and books into a State snapshot.  Note: this
         * acquires a run lock on the simulation, and so will block if the simulation is currently
         * running (i.e. in another thread).
         */
        State(const std::shared_ptr<eris::Simulation> &sim);

        /// Creates an empty state with t=0, no readers, and no books.
        State() = default;

        /// The simulation period represented by this state
        eris::eris_time_t t = 0;

        /** The simulation boundary.  This is inferred from the boundary of the first reader found
         * in the simulation, the first book if there are no readers, and otherwise remains at its
         * default, NaN.
         */
        double boundary = std::numeric_limits<double>::quiet_NaN();

        /** The simulation dimensions.  This is inferred from the first reader found in the
         * simulation, the first book if there are no readers, and otherwise the default, 0.
         */
        uint32_t dimensions = 0;

        /// The readers at the given state
        std::map<eris::eris_id_t, ReaderState> readers;

        /// The books at the given state
        std::map<eris::eris_id_t, BookState> books;

        /// The public tracker, if it exists.
        std::unique_ptr<PublicTrackerState> publictracker;
};

}}
