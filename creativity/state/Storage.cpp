#include "creativity/state/Storage.hpp"

namespace creativity { namespace state {

bool Storage::empty() const { return size() == 0; }
void Storage::reserve(size_t) {}
double Storage::boundary() const { return boundary_; }
uint64_t Storage::sharingBegins() const { return sharing_begins_; }
unsigned int Storage::dimensions() const { return dimensions_; }
void Storage::flush() {}

void Storage::sharingBegins(uint64_t t) {
    if (sharing_begins_ != std::numeric_limits<uint64_t>::max() and t != sharing_begins_)
        throw std::logic_error("Cannot changing storage.sharingBegins()");
    sharing_begins_ = t;
}

}}
