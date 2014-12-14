#pragma once
#include <iterator>
#include <boost/iterator/iterator_facade.hpp>
#include "creativity/state/State.hpp"
#include "creativity/CreativitySettings.hpp"
#include "creativity/state/StorageBackend.hpp"

namespace creativity { namespace state {

/** Base class for state storage which accesses State values.
 *
 * Subclasses must implement the subscript [] operator, size() method, and push_back(&&) method.
 *
 * \sa creativity::state::MemoryStorage
 * \sa creativity::state::FileStorage
 * \sa creativity::state::PsqlStorage
 */
class Storage final {
    public:
        Storage() = delete;

        /** Storage constructor, which can construct a StorageBackend on the fly.
         *
         * \param settings the CreativitySettings reference
         * \param args parameters to forward to the SB constructor
         *
         * \returns a shared pointer to the created Storage object.
         */
        template<class SB, class... Args>
        static typename std::enable_if<std::is_base_of<StorageBackend, SB>::value, std::shared_ptr<Storage>>::type
        create(CreativitySettings &settings, Args&&... args) {
            return std::shared_ptr<Storage>(new Storage(settings, new SB(std::forward<Args>(args)...)));
        }

        /** Returns the State for the given simulation time period.
         *
         * \param t the simulation period desired.
         *
         * \returns A shared pointer to the requested State.  The pointed-at object will be cached,
         * but may only weakly; thus this method may require the backend storage object to load the
         * data from the storage medium on a subsequent call for the same index if the returned
         * reference isn't stored anywhere.
         *
         * \throws std::out_of_range if `t >= size()`
         */
        std::shared_ptr<const State> operator[](eris::eris_time_t t) const;

        /// Returns the number of states currently stored.
        size_t size() const;

        /** Reserves the requested number of states.  By default this does nothing; subclasses
         * should override if they have a useful reserve() implementation.
         */
        void reserve(size_t capacity);

        /// Returns true if the container is empty.
        bool empty() const;

        /** Adds a state to this storage container.  Depending on the backend, the object may only
         * be queued for addition by a thread and not immediately written to the storage medium.
         *
         * \sa flush()
         */
        void push_back(std::shared_ptr<const State> s);

        /** Writes the settings stored in `.settings` (which is a reference to the
         * CreativitySettings given during construction) to the storage medium (if appropriate).
         * Existing settings are replaced.
         *
         * If not called explicitly, this method will be called by the first push_back() call.
         */
        void updateSettings();

        /// Constructs a State using the given arguments, wraps it in a shared_ptr, then inserts it by calling push_back
        template<class... Args>
        void emplace_back(Args&&... args);

        /** Flushes changes of the backend storage object.  This typically blocks until all data in
         * the queue has been written to the underlying storage medium.
         *
         * flush() is called automatically during object destruction, unless flush_on_destroy has
         * been set to false.
         */
        void flush();

        /// Accesses the underlying backend storage instance
        StorageBackend& backend();

        /** If true (the default), flush() is called when the object is destroyed.  If false, it is
         * not, which may result in data loss.
         */
        bool flush_on_destroy = true;

        /// Destructor.  Calls flush() (unless `flush_on_destroy` has been set to false).
        ~Storage();

        /// Random access iterator class for iterating through states
        class state_iterator : public boost::iterator_facade<state_iterator, const std::shared_ptr<const State>, boost::random_access_traversal_tag> {
            private:
                const Storage &storage;
                size_t i = 0;
                std::shared_ptr<const State> curr;

                friend class Storage;
                friend class boost::iterator_core_access;

                state_iterator(const Storage &st, size_t at);

                reference dereference() const;
                void increment();
                void decrement();
                void advance(difference_type n);
                difference_type distance_to(const state_iterator &it) const;
                bool equal(const state_iterator &other) const;
        };

        /** Returns a Random Access Iterator to the beginning of the storage object's states. */
        state_iterator begin() const;
        /** Returns a Random Access Iterator to the just-past-the-end of the storage object's
         * states. */
        state_iterator end() const;

    private:
        /** Storage constructor, which must be called with a CreativitySettings reference and
         * StorageBackend object.
         *
         * \param settings the CreativitySettings reference associated with the Creativity object
         * this storage class represents.
         * \param backend a StorageBackend-derived object
         */
        Storage(CreativitySettings &settings, StorageBackend *sb) : settings_(settings), backend_(sb)
        {
            backend_->readSettings(settings_);
            // Track the number of states using the size of the cache_
            cache_.reserve(backend_->size());
        }

        /** The simulation settings reference, which subclasses must set during construction, typically from a with default initialization of fields to 0.  These settings
         * should be updated as soon as possible.
         */
        CreativitySettings &settings_;

        /** Weak pointers to stored states.  This cache is automatically managed by this base class.
         */
        mutable std::vector<std::weak_ptr<const State>> cache_;

        /** The storage backend, providing the actual low-level storage access. */
        std::unique_ptr<StorageBackend> backend_;

        /** Tracks whether updateSettings() has been called, so that it can be called during
         * push_back (if needed).
         */
        bool need_settings_updated_ = true;
};

template <class... Args>
void Storage::emplace_back(Args&&... args) {
    push_back(std::make_shared<State>(std::forward<Args>(args)...));
}

}}
