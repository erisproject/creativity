#pragma once
#include <queue>
#include <condition_variable>
#include <cstddef>
#include <eris/types.hpp>
#include <memory>
#include <mutex>
#include <thread>
#include <eris/noncopyable.hpp>
#include "creativity/state/State.hpp"

namespace creativity { struct CreativitySettings; }

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
        /** Constructor: takes a CreativitySettings reference.  The reference must stay valid for
         * the duration of the StorageBackend.
         */
        StorageBackend(CreativitySettings &settings) : settings_{settings} {}

        /** Returns the number of states stored in the storage object, not including queued but
         * not-yet-added states.  This is typically called just once by Storage immediately after
         * constructing the object to determine the number of initial storage states; Storage
         * internally tracks (by the number of enqueue calls) the number of states after that.
         */
        virtual size_t size() = 0;

        /// Returns true if the container is empty; exactly equivalent to `size() == 0`
        bool empty();

        /** Optional method that subclasses can override if knowing the number of states that will
         * be stored in advance is useful.  The default implementation does nothing.
         */
        virtual void reserve(size_t T);

        /** Writes the settings in settings_ to the storage medium (if appropriate).  Existing
         * stored settings should be replaced with the new values.
         */
        virtual void writeSettings() = 0;

        /** Flushes changes.  This method blocks until the current thread has finished writing all
         * queued states and the device has been flushed to storage (if applicable and the given
         * argument is true).
         *
         * If a subclass supports additional flushing capabilities (e.g. syncing a file to disk), it
         * should override storage_flush(), which this method calls.
         */
        void flush(bool flush_buffers = true);

        /** Attempts to flush changes, waiting at most the given number of milliseconds for the
         * flush to complete.  Returns true if the flush completed, false if the duration expired
         * before the flush completed.
         *
         * \returns true if the flush completed, false if the duration expired without the flush
         * finishing.
         */
        bool flush_for(long milliseconds, bool flush_buffers = true);

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
        virtual std::shared_ptr<const State> load(eris::eris_time_t t) = 0;

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

        /** Flush any buffered data to the underlying storage system.  The default does nothing;
         * storage backends with such a concept should override.
         */
        virtual void storage_flush();

        /** CreativitySettings reference. */
        CreativitySettings &settings_;

    private:

        /** Queued values awaiting storage in the storage medium.  If a state is not found in
         * cache_, this is checked next.
         */
        std::queue<std::shared_ptr<const State>> queue_;

        /** Number of states values that have been removed from queue_ but are still being written.
         */
        unsigned queue_pending_ = 0;

        /** Set to true in the parent to tell the thread to flush; set to false in the thread once
         * the flush is completed.
         */
        bool flush_pending_ = false;

        /** Mutex guarding access to queue_, queue_pending_, and flush_pending_.  If the storage
         * object uses a thread, it must lock this mutex before manipulating these variables (the
         * methods of this base class automatically obtain such a lock as needed).
         */
        mutable std::mutex queue_mutex_;

        /** Condition variable for the queue for use by a subclass.  This is signalled after the
         * base class adds to queue_.
         */
        std::condition_variable queue_cv_;

        /** Condition variable for the thread to signal that it has emptied the queue and (possibly)
         * flushed the device. Used during flush(). */
        std::condition_variable queue_finished_cv_;

        /** Thread object (for subclasses that return true for threaded()).  The thread runs perpetually
         * until the StorageBackend destructor is called.
         */
        std::thread thread_;

        /** Set during destruction so that the thread loop ends.  Note that if flushing on
         * destruction is desired, flush() on the backend should be called explicitly during the
         * destruction of the object that owns the StorageBackend (i.e. Storage).  The thread
         * behaviour is to abort as soon as this is seen, even if other states/flush requests are
         * pending. */
        bool thread_quit_ = false;

        /** Thread loop.  Calls thread_insert() as needed. */
        void thread_loop_();

};

}}
