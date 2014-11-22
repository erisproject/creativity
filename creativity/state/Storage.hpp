#pragma once
#include <iterator>
#include "creativity/state/State.hpp"
#include "creativity/CreativitySettings.hpp"

namespace creativity { namespace state {

/** Base class for state storage which accesses State values.
 *
 * Subclasses must implement the subscript [] operator, size() method, and push_back(&&) method.
 *
 * \sa creativity::state::MemoryStorage
 * \sa creativity::state::FileStorage
 */
class Storage {
    public:
        /** Returns the State at the given position.  It is highly recommended that subclasses store
         * objects in simulation order, that is, that `storage[j].t == j` is true.
         *
         * \param i the index (if the subclass follows the recommendation above, also the simulation
         * period).
         *
         * \returns A shared pointer to the requested State.  The pointed-at object may or may not
         * be stored internally by the storage object; in particular, storage-based implementations
         * may construct and return a new State object each time this operator is called.
         *
         * \throws std::out_of_range if `t >= size()`
         */
        virtual std::shared_ptr<const State> operator[](size_t i) const = 0;

        /** Accesses the simulation settings associated with this storage object.
         */
        const CreativitySettings &settings = settings_;

        /// Returns the number of states currently stored.
        virtual size_t size() const = 0;

        /** Reserves the requested number of states.  By default this does nothing; subclasses
         * should override if they have a useful reserve() implementation.
         */
        virtual void reserve(size_t capacity);

        /// Returns true if the container is empty.
        virtual bool empty() const;

        /** Adds a state to this storage container.  Note that it is up to the implementing class
         * whether it stores the given std::shared_ptr or just uses it to access the State.
         */
        virtual void push_back(std::shared_ptr<const State> s) = 0;

        /// Constructs a State using the given arguments, wraps it in a shared_ptr, then inserts it by calling push_back
        template<class... Args>
        void emplace_back(Args&&... args);

        /** Flushes changes, if the underlying storage object has such a concept.  The default
         * implementation does nothing.  If this is not called at the end of a program, written data
         * may not actually be saved to the underlying storage system.
         */
        virtual void flush();

        /// Default destructor
        virtual ~Storage() = default;

    protected:
        /** The simulation settings, with default initialization of fields to 0.  These settings
         * should be updated as soon as possible.
         */
        CreativitySettings settings_ = {};
};

template <class... Args>
void Storage::emplace_back(Args&&... args) {
    push_back(std::make_shared<State>(std::forward<Args>(args)...));
}

}}
