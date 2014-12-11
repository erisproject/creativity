#include "creativity/state/Storage.hpp"

namespace creativity { namespace state {

bool Storage::empty() const { return size() == 0; }
void Storage::reserve(size_t) {}
void Storage::flush() {}


void Storage::push_back(std::shared_ptr<const State> state) {
    // Make sure the state agrees with the states we already have
    if (settings_.dimensions != 0 and state->dimensions != settings_.dimensions)
        throw std::runtime_error("Cannot add state: state dimensions differ from storage dimensions");
    if (settings_.boundary != 0 and state->boundary != settings_.boundary)
        throw std::runtime_error("Cannot add state: state boundary differs from storage boundary");
    if (settings_.readers != 0 and settings_.readers != state->readers.size())
        throw std::runtime_error("Cannot add state: state #readers differs from storage #readers");

    if (need_settings_updated_) {
        updateSettings();
        need_settings_updated_ = false;
    }

    push_back_(std::move(state));
}

}}
