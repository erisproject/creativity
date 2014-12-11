#pragma once
#include <iterator>
#include <boost/iterator/iterator_facade.hpp>
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
        Storage() = delete;

        /** Storage base class constructor, which must be called with a CreativitySettings
         * reference.
         */
        Storage(CreativitySettings &settings) : settings_(settings) {}

        /** Returns the State at the given position.  It is highly recommended that subclasses store
         * objects in simulation order, that is, that `storage[j]->t == j` is true.
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
         *
         * This calls the subclass-specific push_back_() method.
         */
        void push_back(std::shared_ptr<const State> s);

        /** Writes the settings stored in `.settings` (which is a reference to the
         * CreativitySettings given during construction) to the storage medium (if appropriate).
         * Existing settings are replaced.
         *
         * Subclasses should ensure that this is called during the first push_back() call if it
         * hasn't yet been called and a new data entry was created (rather than an existing one
         * loaded).
         */
        virtual void updateSettings() = 0;

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

        /// Random access iterator class for iterating through states
        class state_iterator : public boost::iterator_facade<state_iterator, const std::shared_ptr<const State>, boost::random_access_traversal_tag> {
            private:
                const Storage &storage;
                size_t i = 0;
                std::shared_ptr<const State> curr;

                friend class Storage;
                friend class boost::iterator_core_access;

                state_iterator(const Storage &st, size_t at) : storage(st) { advance(at); }

                reference dereference() const { return curr; }
                void increment() { advance(1); }
                void decrement() { advance(-1); }
                void advance(difference_type n) { i += n; if (i < storage.size()) curr = storage[i]; else curr.reset(); }
                difference_type distance_to(const state_iterator &it) { return it.i - i; }
                bool equal(const state_iterator &other) const { return other.i == i; }
        };

        /** Returns a Random Access Iterator to the beginning of the storage object's states. */
        state_iterator begin() const {
            return state_iterator(*this, 0);
        }
        /** Returns a Random Access Iterator to the just-past-the-end of the storage object's
         * states. */
        state_iterator end() const {
            return state_iterator(*this, size());
        }

    protected:
        /** The simulation settings reference, which subclasses must set during construction, typically from a with default initialization of fields to 0.  These settings
         * should be updated as soon as possible.
         */
        CreativitySettings &settings_;

        /** Adds a state to this storage container.  Note that it is up to the implementing class
         * whether it stores the given std::shared_ptr or just uses it to access the State.
         */
        virtual void push_back_(std::shared_ptr<const State> &&s) = 0;

        /** True if only default settings have been written, in which case updateSettings() will be
         * called during push_back() before invoking the subclassable push_back_() method.
         */
        bool need_settings_updated_ = false;
};

template <class... Args>
void Storage::emplace_back(Args&&... args) {
    push_back(std::make_shared<State>(std::forward<Args>(args)...));
}

}}
