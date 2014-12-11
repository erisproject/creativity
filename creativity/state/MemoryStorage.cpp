#include "creativity/state/MemoryStorage.hpp"

namespace creativity { namespace state {

MemoryStorage::MemoryStorage(CreativitySettings &set) : Storage(set) {}

MemoryStorage::MemoryStorage(const Storage &copy) : Storage(copy) {
    reserve(copy.size());
    for (size_t i = 0; i < copy.size(); i++) {
        states_.push_back(copy[i]);
    }
}

std::shared_ptr<const State> MemoryStorage::operator[](size_t i) const { return states_[i]; }

size_t MemoryStorage::size() const { return states_.size(); }

void MemoryStorage::reserve(size_t capacity) { states_.reserve(capacity); }

void MemoryStorage::push_back_(std::shared_ptr<const State> &&s) {
    states_.push_back(std::move(s));
}

}}
