#include "creativity/state/MemoryStorage.hpp"

namespace creativity { namespace state {

MemoryStorage::MemoryStorage(const CreativitySettings &set) {
    settings_ = set;
}

MemoryStorage::MemoryStorage(const Storage &copy) {
    settings_ = copy.settings;
    reserve(copy.size());
    for (size_t i = 0; i < copy.size(); i++) {
        states_.push_back(copy[i]);
    }
}

std::shared_ptr<const State> MemoryStorage::operator[](size_t i) const { return states_[i]; }

size_t MemoryStorage::size() const { return states_.size(); }

void MemoryStorage::reserve(size_t capacity) { states_.reserve(capacity); }

void MemoryStorage::push_back(std::shared_ptr<const State> s) {
    if (settings.boundary != s->boundary) throw std::runtime_error("Cannot add State with different boundaries!");
    if (settings.dimensions != s->dimensions) throw std::runtime_error("Cannot add State with different dimensions!");
    if (settings.readers != s->readers.size()) throw std::runtime_error("Cannot add State with different #readers!");

    states_.push_back(std::move(s));
}

}}
