#include "creativity/state/Storage.hpp"

namespace creativity { namespace state {

bool Storage::empty() const { return size() == 0; }
void Storage::reserve(size_t) {}
double Storage::boundary() const { return boundary_; }
unsigned int Storage::dimensions() const { return dimensions_; }

}}
