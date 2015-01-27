#pragma once
#include <queue>
#include <iterator>
#include <boost/iterator/iterator_facade.hpp>
#include <eris/noncopyable.hpp>
#include "creativity/state/State.hpp"
#include "creativity/CreativitySettings.hpp"

namespace creativity { namespace state {

/** Base class for state storage which accesses State values.
 *
 * Subclasses must set size_ when initialized; after initialization, this class will take care of
 * updating it as new states are added.
 *
 * Subclasses must implement the subscript [] operator, size() method, and push_back(&&) method.
 *
 * \sa creativity::state::MemoryStorage
 * \sa creativity::state::FileStorage
 */
class StorageBackend : private eris::noncopyable {
    public:
        StorageBackend() = default;

        /** Returns the number of states stored in the storage object, not including queued but
         * not-yet-added states.  This is typically called just once by Storage immediately after
         * constructing the object to determine the number of initial storage states; Storage
         * internally tracks (by the number of enqueue calls) the number of states after that.
         */
        virtual size_t size() const = 0;

        /// Returns true if the container is empty, i.e. if `size() == 0`
        bool empty() const;

        /** Subclasses should set this to true if they have some settings to load.  This value is
         * checked immediately after construction by Storage, and should return true if the backend
         * has existing settings to load, false if the backend represented a new creativity storage
         * target without meaningful settings.
         *
         * The default value is false.
         */
        bool have_settings = false;

        /** Optional method that subclasses can override if knowing the number of states that will
         * be stored in advance is useful.  The default implementation does nothing.
         */
        virtual void reserve(size_t T);

        /** Writes the given settings to the storage medium (if appropriate).  Existing settings are
         * replaced.
         */
        virtual void writeSettings(const CreativitySettings &settings) = 0;

        /** Reads the currently stored settings from the storage medium, updating the given object.
         * This is typically called by Storage immediately after backend object construction if
         * hasSettings() returned true.  If no settings are available, this method should just leave
         * the given CreativitySettings object untouched.
         */
        virtual void readSettings(CreativitySettings &settings) const = 0;

        /** Flushes changes.  This method blocks until the current thread has finished writing all
         * queued states.  If there is no thread currently active, this returns immediately.
         *
         * If a subclass supports additional flushing capabilities (e.g. syncing a file to disk), it
         * should override, call the base flush method, then perform the additional flush
         * operations.
         */
        virtual void flush();

        /** Adds a state to this storage container.  The default implementation spawns a thread then
         * adds it to queue_.  Subclasses not using threaded storage MUST override this method to
         * perform the storage immediately.
         */
        virtual void enqueue(std::shared_ptr<const State> &&s);

        /** Returns the number of states still pending in the queue. */
        virtual size_t pending() const;

        /** Gets the requested state from the backend.  This is not required to check for objects in
         * the queue: the caller is responsible for caching queued objects.
         *
         * If the state does not exist in the storage medium, a subclass should return a void
         * shared_ptr.
         */
        virtual std::shared_ptr<const State> load(eris::eris_time_t t) const = 0;

        /// Virtual destructor that tells the thread (if created) to quit and waits for it.
        virtual ~StorageBackend();

    protected:
        /** Called from the queue thread to write the given State to the underlying storage.
         * Threading subclasses MUST override this: the base class implementation throws an
         * exception if called.
         *
         * \throws std::logic_error if called: threading subclasses must override this method.
         */
        virtual void thread_insert(std::shared_ptr<const State> &&s);

    private:

#ifndef CREATIVITY_DISABLE_THREADED_STORAGE

        /** Queued values awaiting storage in the storage medium.  If a state is not found in
         * cache_, this is checked next.
         */
        std::queue<std::shared_ptr<const State>> queue_;

        /** Mutex guarding access to queue_.  If the storage object uses a thread, it must lock this
         * mutex before manipulating queue_.  The methods of this base class automatically obtain a
         * lock before using queue_.
         */
        mutable std::mutex queue_mutex_;

        /** Condition variable for the queue for use by a subclass.  This is signalled after the
         * base class adds to queue_.
         */
        std::condition_variable queue_cv_;

        /** Condition variable for the thread to signal that it has emptied the queue. Used during
         * flush(). */
        std::condition_variable queue_finished_cv_;

        /** Thread object (for subclasses that return true for threaded()).  The thread loops perpetually
         * until the StorageBackend destructor is called.
         */
        std::thread thread_;

        /** Set during destruction so that the thread loop ends. */
        bool thread_quit_ = false;

        /** Thread loop.  Calls thread_insert() as needed. */
        void thread_inserter_();
#endif

};

}}
