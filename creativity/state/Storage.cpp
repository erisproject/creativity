#include "creativity/state/Storage.hpp"

using namespace eris;

namespace creativity { namespace state {

size_t Storage::size() const { return cache_.size(); }
bool Storage::empty() const { return cache_.empty(); }
void Storage::reserve(size_t T) { backend_->reserve(T); }
void Storage::flush() { backend_->flush(); }

StorageBackend& Storage::backend() { return *backend_; }

Storage::~Storage() { if (flush_on_destroy) flush(); }

Storage::state_iterator::state_iterator(const Storage &st, size_t at) : storage(&st) { advance(at); }
Storage::state_iterator::reference Storage::state_iterator::dereference() const { return curr; }
void Storage::state_iterator::increment() { advance(1); }
void Storage::state_iterator::decrement() { advance(-1); }
void Storage::state_iterator::advance(difference_type n) {
    i += n;
    if (i < storage->size()) curr = (*storage)[i];
    else curr.reset();
}
Storage::state_iterator::difference_type Storage::state_iterator::distance_to(const state_iterator &it) const {
    return it.i - i;
}
bool Storage::state_iterator::equal(const state_iterator &other) const { return other.i == i; }
Storage::state_iterator Storage::begin() const { return state_iterator(*this, 0); }
Storage::state_iterator Storage::end() const { return state_iterator(*this, size()); }

std::shared_ptr<const State> Storage::operator[](eris_time_t t) const {
    // First try the cache:
    std::shared_ptr<const State> st = cache_[t].lock();
    if (not st) {
        // Not found; try loading from the backend:
        st = backend_->load(t);
        // If that didn't work, it doesn't exist, so throw:
        if (not st) throw std::out_of_range("StorageBackend::operator[]: requested state does not exist in storage object");

        // Cache the value before returning it
        if (cache_.size() <= t) cache_.resize(t+1);
        cache_[t] = st;
    }

    return st;
}

void Storage::push_back(std::shared_ptr<const State> state) {
    // Make sure the state agrees with the states we already have
    if (settings_.dimensions != 0 and state->dimensions != settings_.dimensions)
        throw std::runtime_error("Cannot add state: state dimensions differ from storage dimensions");
    if (settings_.boundary != 0 and state->boundary != settings_.boundary)
        throw std::runtime_error("Cannot add state: state boundary differs from storage boundary");

    if (need_settings_updated_) {
        updateSettings();
    }

    cache_.emplace_back(state);
    backend_->enqueue(std::move(state));
}

void Storage::updateSettings() {
    backend_->writeSettings(settings_);
    need_settings_updated_ = false;
}

}}
