#include "creativity/state/StorageBackend.hpp"
#include <algorithm>
#include <stdexcept>

using namespace eris;

namespace creativity { namespace state {

bool StorageBackend::empty() const { return size() == 0; }

void StorageBackend::reserve(size_t) {}

size_t StorageBackend::pending() const {
    if (not thread_.joinable()) return 0;
    std::unique_lock<std::mutex> lock(queue_mutex_);
    return queue_.size() + queue_pending_;
}

void StorageBackend::enqueue(std::shared_ptr<const State> &&s) {
    if (not thread_.joinable())
        thread_ = std::thread(&StorageBackend::thread_inserter_, this);

    {
        std::unique_lock<std::mutex> lock(queue_mutex_);
        queue_.push(std::move(s));
    }
    queue_cv_.notify_one();
}

void StorageBackend::thread_insert(std::shared_ptr<const State>&&) {
    throw std::logic_error("StorageBackend::thread_queue called: a subclass forgot to override either enqueue or thread_insert!");
}

void StorageBackend::thread_inserter_() {
    std::unique_lock<std::mutex> lock(queue_mutex_);

    while (true) {
        queue_cv_.wait(lock, [&] { return thread_quit_ or not queue_.empty(); });

        while (not thread_quit_ and not queue_.empty()) {
            // Pop off one at a time here so that if a thread_quit_ comes in, we'll respond (without
            // having to wait to write out the entire queue).
            auto st = queue_.front();
            queue_.pop();
            queue_pending_++;

            lock.unlock();
            thread_insert(std::move(st));
            lock.lock();
            queue_pending_--;
        }

        // If the queue is empty now, send a signal in case the main thread is waiting (during
        // a flush)
        if (queue_.empty()) queue_finished_cv_.notify_one();

        // If asked to quit, do so:
        if (thread_quit_) return;
    }
}

void StorageBackend::flush() {
    if (thread_.joinable()) {
        std::unique_lock<std::mutex> lock(queue_mutex_);
        queue_finished_cv_.wait(lock, [&] { return queue_.empty() and queue_pending_ == 0; });
    }
}

bool StorageBackend::flush_for(long milliseconds) {
    if (thread_.joinable()) {
        std::unique_lock<std::mutex> lock(queue_mutex_);
        return queue_finished_cv_.wait_for(lock, std::chrono::milliseconds(milliseconds), [&] { return queue_.empty() and queue_pending_ == 0; });
    }
    return true;
}

StorageBackend::~StorageBackend() {
    if (thread_.joinable()) {
        {
            std::unique_lock<std::mutex> lock(queue_mutex_);
            thread_quit_ = true;
        }
        queue_cv_.notify_one();
        thread_.join();
    }
}
    

}}
