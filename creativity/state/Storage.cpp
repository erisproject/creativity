#include "creativity/state/Storage.hpp"

namespace creativity { namespace state {

bool Storage::empty() const { return size() == 0; }
void Storage::reserve(size_t) {}
void Storage::flush() {}

}}
