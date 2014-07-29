/** \file
 *
 * This file contains the inline read_value, write_value, and parse_value implementations for a
 * big-endian system.  Since the file storage is little-endian, byte orders have to be reversed when
 * converting between byte strings and memory locations.
 *
 * This file is loaded automatically on big-endian systems from FileStorage.hpp.
 */
#pragma once

// These are documented in FileStorage.hpp, so tell doxygen to ignore:
/// \cond
template <typename T>
inline T FileStorage::read_value() const {
    T val;
    char *valbytes = reinterpret_cast<char*>(&val);
    char input[sizeof(T)];
    f_.read(input, sizeof(T));
    for (size_t i = 0; i < sizeof(T); i++) valbytes[i] = input[sizeof(T) - 1 - i];
    return val;
}

template <typename T>
inline void FileStorage::write_value(const T &val) {
    const char *in = reinterpret_cast<const char*>(&val);
    char out[sizeof(T)];
    for (size_t i = 0; i < sizeof(T); i++) out[i] = in[sizeof(T) - 1 - i];
    f_.write(out, sizeof(T));
}

template <typename T, typename F>
inline T FileStorage::parse_value(const F &from) {
    T val;
    const char *from_bytes = reinterpret_cast<const char*>(&from);
    char *valbytes = reinterpret_cast<char*>(&val);
    for (size_t i = 0; i < sizeof(T); i++) valbytes[i] = from_bytes[sizeof(T) - 1 - i];
    return val;
}

template <typename F, typename T>
inline void FileStorage::store_value(const F &from, T &to) {
    const char *from_bytes = reinterpret_cast<const char*>(&from);
    char *to_bytes = reinterpret_cast<char*>(&to);
    for (size_t i = 0; i < sizeof(F); i++) to_bytes[i] = from_bytes[sizeof(F) - 1 - i];
}

/// \endcond

