#include "creativity/state/MemoryStorage.hpp"

using namespace eris;

namespace creativity { namespace state {

MemoryStorage::MemoryStorage(const Storage &copy) {
    reserve(copy.size());
    states_.insert(states_.end(), copy.begin(), copy.end());
}

void MemoryStorage::writeSettings(const CreativitySettings&) {}
void MemoryStorage::readSettings(CreativitySettings&) const {}

std::shared_ptr<const State> MemoryStorage::load(eris_time_t i) const {
    if (i < states_.size()) return states_[i];
    return std::shared_ptr<const State>();
}

size_t MemoryStorage::size() const { return states_.size(); }

void MemoryStorage::reserve(size_t capacity) { states_.reserve(capacity); }

void MemoryStorage::enqueue(std::shared_ptr<const State> &&s) {
    states_.push_back(std::move(s));
}

}}
