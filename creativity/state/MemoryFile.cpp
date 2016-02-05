#include "creativity/state/MemoryFile.hpp"

namespace creativity { namespace state {

MemoryFile::MemoryFile(const std::string &data, MODE mode) {
    std::string blank;
    const std::string &initial_data = mode == MODE::OVERWRITE ? blank : data;
    f_.reset(new std::stringstream(initial_data, mode == MODE::READONLY ? open_readonly : open_readwrite));
    f_->exceptions(f_->failbit | f_->badbit);

    bool empty = (mode == MODE::OVERWRITE ? true :
            (mode == MODE::READ or mode == MODE::READONLY) ? false :
            initial_data.empty());

    initialize(empty);
}

void MemoryFile::storage_flush() {}

}}
