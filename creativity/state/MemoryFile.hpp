#pragma once
#include "creativity/state/FileStorage.hpp"
#include <sstream>

namespace creativity { namespace state {

/** Class that extends FileStorage to write to an in-memory buffer instead of a file.  This means
 * that the results will vanish as soon as this object vanishes, and that an amount of memory
 * greater (depending on the implementation, but often the next power of 2) the would-be file is
 * used, depending on the stl c++ library implementation of sstream.
 *
 * \sa FileStorage
 */
class MemoryFile final : public FileStorage {
    public:
        /** Creates an in-memory file consisting of a copy of the given string.
         *
         * The mode argument is exactly like those supplied to the FileStorage constructor.  Note,
         * in particular, that this means MODE::OVERWRITE simply ignores the given string.
         *
         * \param data the file data to copy into memory; defaults to an empty string.
         * \param mode the "file" mode; defaults to MODE::APPEND
         * \sa FileStorage()
         */
        explicit MemoryFile(const std::string &data = "", MODE mode = MODE::APPEND);

        /** Does nothing. */
        virtual void storage_flush() override;
};

}}
