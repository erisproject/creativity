/** \file
 *
 * This file contains the inline read_value, write_value, and parse_value implementations for a
 * little-endian system (such as x86/AMD64).  Since the file storage is little-endian, no special
 * conversion is needed: this code simply copies the bytes into the appropriate types, or interprets
 * the type as bytes.
 *
 * This file is loaded automatically on little-endian systems from FileStorage.hpp.
 */
#pragma once

// These are documented in FileStorage.hpp, so tell doxygen to ignore:
/// \cond
template <typename T>
inline T FileStorage::read_value() const {
    T val;
    f_.read(reinterpret_cast<char*>(&val), sizeof val);
    return val;
}

template <typename T>
inline void FileStorage::write_value(const T &val) {
    f_.write(reinterpret_cast<const char*>(&val), sizeof(T));
}

template <typename T, typename F>
inline T FileStorage::parse_value(const F &from) {
    T val;
    // If F and T are the same type (and thus have the same alignment requirements) we don't have to
    // do any extra work and can simply copy by value.  This is only valid when storage and memory
    // orders are the same, however, and so other endian storage methods can't do this.
    if (std::is_same<T, F>::value) {
        val = from;
    }
    // Otherwise we have to copy the bytes starting at the given value into a new variable of the
    // requested type, to ensure that we don't have alignment issues (for example, `from` could be a
    // char with memory alignment `(&from % 8 == 5)`, but `T` could be a double requiring 8-byte
    // alignment, and so interpreting the memory starting at &from as a T won't be valid.
    else {
        const char *from_bytes = reinterpret_cast<const char*>(&from);
        std::copy(&from_bytes[0], &from_bytes[sizeof(T)], reinterpret_cast<char*>(&val));
    }
    return val;
}

template <typename F, typename T>
inline void FileStorage::store_value(const F &from, T &to) {
    // If T and F are the same, simply copy-by-value (see comment in parse_value).
    if (std::is_same<T, F>::value) {
        to = from;
    }
    // Otherwise copy the bytes.
    else {
        const char *from_bytes = reinterpret_cast<const char*>(&from);
        std::copy(&from_bytes[0], &from_bytes[sizeof(F)], reinterpret_cast<char*>(&to));
    }
}
    
/// \endcond

