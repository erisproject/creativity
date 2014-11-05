#include "creativity/state/FileStorage.hpp"
#include <cmath>
#include <sys/stat.h>

namespace creativity { namespace state {

using namespace creativity::belief;
using namespace Eigen;
using eris::eris_id_t;

/// Member definitions for header/cblock constants
constexpr unsigned int
        FileStorage::HEADER::size,
        FileStorage::HEADER::states,
        FileStorage::CBLOCK::size,
        FileStorage::CBLOCK::states;
constexpr char FileStorage::HEADER::fileid[], FileStorage::HEADER::test_value[];
constexpr uint64_t FileStorage::HEADER::u64_test;
constexpr int64_t
        FileStorage::HEADER::pos::fileid,
        FileStorage::HEADER::pos::filever,
        FileStorage::HEADER::pos::test_value,
        FileStorage::HEADER::pos::states,
        FileStorage::HEADER::pos::dimensions,
        FileStorage::HEADER::pos::boundary,
        FileStorage::HEADER::pos::state_first,
        FileStorage::HEADER::pos::state_last,
        FileStorage::HEADER::pos::continuation,
        FileStorage::HEADER::i64_test,
        FileStorage::CBLOCK::next_cblock;
constexpr uint32_t FileStorage::HEADER::u32_test;
constexpr int32_t FileStorage::HEADER::i32_test;
constexpr double FileStorage::HEADER::dbl_test;
constexpr std::ios_base::openmode FileStorage::open_readonly, FileStorage::open_overwrite, FileStorage::open_readwrite;

FileStorage::ParseError::ParseError(const std::string &message) : std::runtime_error(message) {}

void FileStorage::throwParseError(const std::string& message) const {
    decltype(f_.tellg()) pos;
    try { pos = f_.tellg(); }
    catch (std::ios_base::failure &f) { pos = -1; } // -1 is also be returned by .tellg() for pre-existing failures
    throw ParseError("Parsing file failed [pos=" + (pos == -1 ? std::string{"(error)"} : std::to_string(pos)) + "]: " + message);
}

inline bool file_exists(const std::string &name) {
    struct stat buffer;
    return stat(name.c_str(), &buffer) == 0;
}

FileStorage::FileStorage(const std::string &filename, MODE mode) {
    f_.exceptions(f_.failbit | f_.badbit);

    bool parse = true;
    switch (mode) {
        case MODE::READONLY:
            f_.open(filename, open_readonly);
            parse = true;
            break;
        case MODE::READ:
            f_.open(filename, open_readwrite);
            parse = true;
            break;
        case MODE::OVERWRITE:
            f_.open(filename, open_overwrite);
            parse = false;
            break;
        case MODE::APPEND:
            if (file_exists(filename)) {
                f_.open(filename, open_readwrite);
                f_.seekg(0, f_.end);
                parse = f_.tellp() > 0; // Parse if the file is non-empty, write otherwise
            }
            else {
                // Open if overwrite mode if the file doesn't exist (open in read-write mode
                // fails if the file doesn't exist (unless also truncating)).
                f_.open(filename, open_overwrite);
                parse = false;
            }
            break;
    }

    if (parse)
        parseMetadata();
    else
        writeEmptyHeader();
}

size_t FileStorage::size() const {
    return state_pos_.size();
}

void FileStorage::push_back(std::shared_ptr<const State> state) {
    // Make sure the state agrees with the states we already have
    if (dimensions_ != 0 and state->dimensions != dimensions_)
        throw std::runtime_error("Cannot add state: state dimensions differ from storage dimensions");
    if (boundary_ != 0 and state->boundary != boundary_)
        throw std::runtime_error("Cannot add state: state boundary differs from storage boundary");

    f_.seekp(0, f_.end);
    auto location = f_.tellp();
    if (location < HEADER::size) {
        throw std::runtime_error("File size is invalid, unable to write new State");
    }
    // We're now at the end of the file, which is where we want to put it.

    // Write out the record before storing its location: if writing fails, there may be some
    // unreferences garbage (a partially written state) at the end of the file, but the file will
    // still be consistent: the garbage won't be referenced anywhere and so will just be ignored
    // if the file is reloaded.
    writeState(*state);

    addStateLocation(location);

    if (states_.size() < state_pos_.size()) states_.resize(state_pos_.size());
    states_[state_pos_.size()-1] = state;

    if (dimensions_ == 0 and state->dimensions != 0) {
        f_.seekp(HEADER::pos::dimensions);
        write_u32(state->dimensions);
        dimensions_ = state->dimensions;
    }
    if (boundary_ == 0 and state->boundary != 0) {
        f_.seekp(HEADER::pos::boundary);
        write_value(state->boundary);
        boundary_ = state->boundary;
    }
}


std::shared_ptr<const State> FileStorage::operator[](size_t i) const {
    std::shared_ptr<const State> ret;
    if (states_.size() > i and (ret = states_[i].lock())) {
        return ret;
    }

    if (i >= state_pos_.size())
        throw std::out_of_range("state::FileStorage: requested State index is invalid");

    auto state_pos = state_pos_[i];
    f_.seekg(state_pos);

    ret = readState();
    if (states_.size() <= i) states_.resize(i+1);
    states_[i] = ret;
    return ret;
}

void FileStorage::flush() {
    f_.flush();
}

void FileStorage::writeEmptyHeader() {
    f_.seekp(0);
    f_.write(HEADER::fileid, sizeof HEADER::fileid); // 'CrSt' file signature
    write_u32(1); // file version
    f_.write(HEADER::test_value, sizeof HEADER::test_value); // All test values squished into one
    write_u32(0); // Number of states (will be updated once known)
    write_u32(0); // Number of dimensions (updated once known)
    write_dbl(0.0); // Boundary (updated once known)

    // State addresses and continuation location:
    char zeros[8 * (HEADER::states + 1)] = {0};
    f_.write(zeros, sizeof zeros);

    if (f_.tellp() != HEADER::size) {
        // If this exception occurs, something in the above sequence is wrong.
        throw std::runtime_error("Header writing failed: header size != " + std::to_string(HEADER::size) + " bytes");
    }
}

void FileStorage::addStateLocation(std::streampos location) {
    uint32_t curr_states = state_pos_.size();
    int past_header = curr_states - HEADER::states;
    if (past_header >= 0 and past_header % CBLOCK::states == 0) {
        // Either the header or the last continuation block was filled up by the last state, so we
        // need to create a new continuation block
        createContinuationBlock();
    }

    if (past_header < 0) {
        f_.seekp(HEADER::pos::state_first + 8 * (int64_t) curr_states);
    }
    else {
        f_.seekp(cont_pos_.back() + (std::streampos) 8 * (past_header % CBLOCK::states));
    }
    write_value(location);

    f_.seekp(HEADER::pos::states);
    write_u32((curr_states + 1));
    state_pos_.push_back(location);
}

void FileStorage::createContinuationBlock() {
    f_.seekp(0, f_.end);
    std::streampos location = f_.tellp();
    char blank[CBLOCK::size] = {0};
    f_.write(blank, CBLOCK::size);
    // Write the new block location either in the header (if this is the first cblock) or in the
    // previous cblock
    f_.seekp(cont_pos_.empty() ? HEADER::pos::continuation : (int64_t) cont_pos_.back() + CBLOCK::next_cblock);
    write_i64(location);
    cont_pos_.push_back(location);
}

void FileStorage::parseMetadata() {
    f_.seekg(0, f_.end);
    auto file_size = f_.tellg();
    if (file_size == 0)
        throwParseError("file contains no data");
    else if (file_size < HEADER::size)
        throwParseError("file is too small to contain header data");

    char block[HEADER::size];
    f_.seekg(0);
    f_.read(block, HEADER::size);

    if (block[0] != 'C' or block[1] != 'r' or block[2] != 'S' or block[3] != 't')
        throwParseError("'CrSt' file signature not found");

    auto version = parse_value<uint32_t>(block[4]);
    if (version != 1)
        throwParseError("encountered unknown version " + std::to_string(version) + " != 1");

    // Okay, so we've got a CrSt version 1 file.  The version number lets the file format change
    // later, if necessary, in which case this code will need to handle the different versions.  For
    // now there is one and only one version, 1.

    // Sanity check: make sure these types are the size we expect
#define CHECK_SIZE(TYPE, SIZE) if (sizeof(TYPE) != SIZE) throwParseError("sizeof(" #TYPE ") != " #SIZE)
    CHECK_SIZE(uint32_t, 4);
    CHECK_SIZE(int32_t, 4);
    CHECK_SIZE(uint64_t, 8);
    CHECK_SIZE(int64_t, 8);
    CHECK_SIZE(std::streamoff, 8);
    CHECK_SIZE(double, 8);
#undef CHECK_SIZE

    // Read the test values and make sure what we interpret equals what should be there (guarding
    // against a potential endian failure).
#define CHECK_TEST(T, TYPE) { \
    auto test = parse_value<TYPE>(block[HEADER::pos::test_value]); \
    if (test != HEADER::T##_test) throwParseError("encoded " #T " test value decoded incorrectly (expected " + \
            std::to_string(HEADER::T##_test) + ", got " + std::to_string(test) + ")"); }
    CHECK_TEST(u32, uint32_t);
    CHECK_TEST(i32, int32_t);
    CHECK_TEST(u64, uint64_t);
    CHECK_TEST(i64, int64_t);
    CHECK_TEST(dbl, double);

    // Now the number of states and number of dimensions:
    auto num_states = parse_value<uint32_t>(block[HEADER::pos::states]);
    dimensions_ = parse_value<uint32_t>(block[HEADER::pos::dimensions]);
    if (dimensions_ == 0 and num_states != 0)
        throwParseError("found invalid dimensions == 0 when num_states > 0");

    boundary_ = parse_value<double>(block[HEADER::pos::boundary]);
    if (boundary_ < 0 or (boundary_ == 0 and num_states != 0))
        throwParseError("found invalid boundary position");

    state_pos_.reserve(num_states);

    size_t header_states = std::min(num_states, HEADER::states);
    parseStateLocations(block[HEADER::pos::state_first], header_states, file_size);

    uint32_t remaining = num_states - header_states;
    while (remaining > 0) {
        auto cont = parse_value<std::streampos>(block[HEADER::pos::continuation]);
        if (cont < HEADER::size or cont >= file_size)
            throwParseError("found invalid continuation location");
        cont_pos_.push_back(cont);

        char cblock[CBLOCK::size];
        f_.seekg(cont, f_.beg);
        f_.read(cblock, sizeof cblock);

        uint32_t block_locations = std::min(remaining, CBLOCK::states);
        parseStateLocations(cblock[0], block_locations, file_size);

        remaining -= block_locations;
    }

    if (num_states != state_pos_.size())
        throwParseError("found " + std::to_string(state_pos_.size()) +
                " state data locations but expected " + std::to_string(num_states));

    // Done!
}

void FileStorage::parseStateLocations(const char &from, const size_t count, const size_t end) {
    const char *data = &from;
    size_t to = count * 8;
    for (size_t i = 0; i < to; i += 8) {
        auto loc = parse_value<std::streamoff>(data[i]);
        if (loc < HEADER::size or (size_t) loc >= end)
            throwParseError("found invalid state data location in header: " + std::to_string(loc));
        state_pos_.push_back(loc);
    }
}

std::shared_ptr<const State> FileStorage::readState() const {
    auto shst = std::make_shared<const State>();
    State &state = const_cast<State&>(*shst);

    state.t = read_u64();
    auto num_readers = read_u32();
    auto num_books = read_u32();

    state.readers.reserve(num_readers);
    state.books.reserve(num_books);

    for (uint32_t i = 0; i < num_readers; i++) {
        state.readers.insert(readReader());
    }

    for (uint32_t i = 0; i < num_books; i++) {
        state.books.insert(readBook());
    }

    return shst;
}

std::pair<eris::eris_id_t, ReaderState> FileStorage::readReader() const {
    auto pair = std::make_pair<eris_id_t, ReaderState>(
            read_u64(),
            ReaderState(dimensions_));

    ReaderState &r = pair.second;
    r.id = pair.first;

    // Position
    for (uint32_t d = 0; d < dimensions_; d++) {
        r.position[d] = read_dbl();
    }
    // Friends
    auto num_friends = read_u32();
    r.friends.reserve(num_friends);
    for (uint32_t i = 0; i < num_friends; i++) {
        r.friends.insert(read_u64());
    }

    // Library
    auto libsize = read_u32();
    r.library.reserve(libsize);
    for (uint32_t i = 0; i < libsize; i++) {
        r.library.emplace(read_u64(), read_dbl());
    }
    // New books
    auto newsize = read_u32();
    r.new_books.reserve(newsize);
    for (uint32_t i = 0; i < newsize; i++) {
        r.new_books.insert(read_u64());
    }
    // Authored books
    auto wrotesize = read_u32();
    r.wrote.reserve(wrotesize);
    for (uint32_t i = 0; i < wrotesize; i++) {
        r.wrote.push_back(read_u64());
    }
    // Utility
    r.u = read_dbl();
    r.u_lifetime = read_dbl();
    // Costs:
    r.cost_fixed = read_dbl();
    r.cost_unit = read_dbl();
    r.income = read_dbl();

    // Beliefs
    belief_data belief = readBelief();
    if (belief.K > 0)
        r.profit = belief.noninformative
            ? Profit(dimensions_, belief.K)
            : Profit(dimensions_, belief.beta, belief.s2, belief.V, belief.n);

    belief = readBelief();
    if (belief.K > 0)
        r.profit_extrap = belief.noninformative
            ? Profit(dimensions_, belief.K)
            : Profit(dimensions_, belief.beta, belief.s2, belief.V, belief.n);

    belief = readBelief();
    if (belief.K > 0)
        r.demand = belief.noninformative
            ? Demand(dimensions_, belief.K)
            : Demand(dimensions_, belief.beta, belief.s2, belief.V, belief.n);

    belief = readBelief();
    if (belief.K > 0)
        r.quality = belief.noninformative
            ? Quality(belief.K)
            : Quality(belief.beta, belief.s2, belief.V, belief.n);

    auto pstream_locs = read_u32();
    for (uint32_t i = 0; i < pstream_locs; i++) {
        belief = readBelief();
        if (belief.K == 0) throwParseError("found illegal profitStream belief with K = 0 (i.e. default constructed model)");
        else if (r.profit_stream.count(belief.K)) throwParseError("found duplicate K value in profit_stream beliefs");
        r.profit_stream.emplace(belief.beta.rows(), belief.noninformative
                ? ProfitStream(belief.K)
                : ProfitStream(belief.beta, belief.s2, belief.V, belief.n));
    }

    return pair;
}

FileStorage::belief_data FileStorage::readBelief() const {
    auto location = read_i64();
    auto return_loc = f_.tellg();
    bool seek_back = false;
    belief_data belief;
    /// Handle magic values:
    if (location == 0) {
        belief.K = 0;
        return belief;
    }
    else if (location < 0 and location >= -100) {
        belief.K = (uint32_t) -location;
        belief.noninformative = true;
        return belief;
    }
    else if (location == -512) {
        // Magic value for "immediately follows".  We need a special value because
        // putting the location there instead means "seek here, read, then seek back"
        // but -512 means "stay here, keep reading, don't seek back" (that is, the encompassing
        // record continues immediately after the belief).

        // Update location to the current location, so that it's the right value when caching
        location = f_.tellg();
    }
    else if (location > HEADER::size) {
        seek_back = true;
        f_.seekg(location);
    }
    else {
        throwParseError("found invalid belief location value");
    }

    auto k = read_u32();
    if (k == 0) throwParseError("Found invalid belief record with k = 0");

    auto size = beliefRecordSize(k);
    // Let std::vector take care of the memory allocation; we'll allocate it with a double value
    // type so that things are aligned properly (which saves parse_value a bit of work, at least on
    // little-endian systems).
    std::vector<double> record(size);
    char *data = reinterpret_cast<char*>(record.data());
    uint32_t *u32data = reinterpret_cast<uint32_t*>(data);

    f_.read(data, size*sizeof(double));

    // Update the belief location cache with the read belief (to enable belief deduplication
    // compression during belief writing).

    // Generate the checksum by treating the little-endian data as native endian u32s.  Since the
    // checksum is only ever stored in memory, it doesn't matter if the file and system endianness
    // differ, so long as checksum calculations are done the same way everywhere.
    uint32_t checksum = std::accumulate(&u32data[0], &u32data[size*(sizeof(double)/sizeof(uint32_t))], uint32_t{0});

    bool already_known = false;
    auto range = belief_locations_[k].equal_range(checksum);
    for (auto it = range.first; it != range.second; ++it) {
        if (it->second == location) {
            already_known = true;
            break;
        }
    }

    if (not already_known) {
        auto val = std::pair<uint32_t, int64_t>(checksum, location);
        if (range.first != belief_locations_[k].end())
            // If our range lookup from before found something, pass the hint along for the new element
            belief_locations_[k].emplace_hint(range.first, std::move(val));
        else
            belief_locations_[k].emplace(std::move(val));
    }

    // The cache record is set, now extract the values from the record block.

    belief.K = k;
    belief.noninformative = false;
    belief.beta = VectorXd(k);
    belief.s2 = 0.0;
    belief.n = 0.0;
    belief.V = MatrixXd(k, k);

    // First K elements are beta values
    uint32_t offset = 0;
    for (; offset < k; ++offset)
        belief.beta[offset] = parse_value<double>(record[offset]);

    // Then s2 and n:
    belief.s2 = parse_value<double>(record[offset++]);
    belief.n = parse_value<double>(record[offset++]);

    // The last K*(K+1)/2 are the V values
    for (unsigned int r = 0; r < k; r++) {
        for (unsigned int c = r; c < k; c++) {
            double cov = parse_value<double>(record[offset++]);
            belief.V(r,c) = cov;
            belief.V(c,r) = cov;
        }
    }

    if (seek_back) f_.seekg(return_loc);

    return belief;
}

void FileStorage::writeState(const State &state) {
    write_u64(state.t);
    write_u32(state.readers.size());
    write_u32(state.books.size());
    for (auto &r : state.readers) {
        writeReader(r.second);
    }
    for (auto &b : state.books) {
        writeBook(b.second);
    }
}

void FileStorage::writeReader(const ReaderState &r) {
    write_u64(r.id);
    for (size_t i = 0; i < r.position.dimensions; i++)
        write_dbl(r.position[i]);

    write_u32(r.library.size());
    for (auto &l : r.library) {
        write_u64(l.first);
        write_dbl(l.second);
    }
    write_u32(r.new_books.size());
    for (auto &n : r.new_books) write_u64(n);
    write_u32(r.wrote.size());
    for (auto &w : r.wrote) write_u64(w);
    write_value(r.u);
    write_value(r.u_lifetime);

    write_value(r.cost_fixed);
    write_value(r.cost_unit);
    write_value(r.income);

    writeBelief(r.profit);
    writeBelief(r.profit_extrap);
    writeBelief(r.demand);
    writeBelief(r.quality);

    // Have to prefix the profit stream beliefs with the count, but since we need to exclude any
    // default-constructed values, we need to loop twice:
    uint32_t ps_count = 0;
    for (auto &b : r.profit_stream) {
        if (b.second.K() > 0) ps_count++;
    }
    write_u32(ps_count);
    for (auto &b : r.profit_stream) {
        if (b.second.K() > 0) writeBelief(b.second);
    }
}

void FileStorage::writeBelief(const Linear &m) {
    auto &k = m.K();
    if (k == 0) {
        // If K = 0 (a default-constructed, non-model object), we just write out a 0 and we're done
        write_i64(0);
        return;
    }
    else if (m.noninformative()) {
        // Non-informative models are easy, too: just write out -K
        write_i64(-(int64_t) k);
        return;
    }

    // Otherwise we need to calculate the data to be written out, then check to see whether it
    // exists in the cache.  If it does, we write the existing location; otherwise we write -512
    // followed by a u32 containing k, then the model values, and store the value (starting at the
    // 'k' value) in the cache.

    // Calculate sizes of the stored data, not including the initial u32 'k' value.
    // `size` is the number of doubles
    auto size = beliefRecordSize(k);
    // `size_u32` is the size when interpreted as u32s
    auto size_u32 = size * (sizeof(double) / sizeof(uint32_t));
    // `size_bytes` is the size when interpreted as chars.
    auto size_bytes = size * sizeof(double);

    // Let std::vector allocate the storage we need.  Note that the values in here might not be the
    // actual double values: the memory will be used to store the values in little-endian order
    // (i.e. the file order), regardless of the system endianness.  In other words, don't read
    // double values out of this vector.
    std::vector<double> record(size);

    size_t pos = 0;
    // First K elements are the beta values
    for (; pos < k; ++pos) {
        store_value(m.beta()[pos], record[pos]);
    }

    // Then s2 and n:
    store_value(m.s2(), record[pos++]);
    store_value(m.n(), record[pos++]);

    // The last k*(k+1)/2 are the lower triangle of the V matrix, in column major order
    for (unsigned int r = 0; r < k; r++) {
        for (unsigned int c = r; c < k; c++) {
            store_value(m.V()(r,c), record[pos++]);
        }
    }

    // Safety check that we wrote the amount we expected to write
    if (pos != size) throw std::runtime_error("Internal error: writeBelief wrote wrong number of elements");

    uint32_t *rec_u32 = reinterpret_cast<uint32_t*>(record.data());
    uint32_t checksum = std::accumulate(&rec_u32[0], &rec_u32[size_u32], uint32_t{0});

    auto curr_loc = f_.tellp();
    bool seek_back = false;
    int64_t found_existing = 0;
    // Allocate space to store the values we need to check:
    std::vector<uint32_t> check_value(size_u32);
    // Now check all matches with the same checksum, read the value, and see if it is exactly
    // equal (it's possible to have duplicate checksums with different data).
    auto range = belief_locations_[k].equal_range(checksum);
    for (auto it = range.first; it != range.second; ++it) {
        auto &loc = it->second;
        seek_back = true;
        f_.seekg(loc + sizeof(uint32_t)); // Skip the "k" value
        f_.read(reinterpret_cast<char*>(check_value.data()), size_bytes);
        if (std::equal(check_value.cbegin(), check_value.cend(), &rec_u32[0])) {
            found_existing = loc;
            break;
        }
    }

    // If we seeked away to check possible matches, seek back:
    if (seek_back) f_.seekp(curr_loc);

    if (found_existing > 0 and false) {
        // If we found an exact match, our record is just its memory location
        write_i64(found_existing);
    }
    else {
        // Otherwise we write the magic value -512 which means "immediately following",
        // then the k, then the record data.  We also store this location in the cache.
        write_i64(-512);
        int64_t cache_loc = f_.tellp();
        write_u32(k);
        f_.write(reinterpret_cast<char*>(record.data()), size_bytes);

        auto val = std::pair<uint32_t, int64_t>(checksum, cache_loc);
        if (range.first != belief_locations_[k].end())
            // If our range found something, pass it along as a hint
            belief_locations_[k].emplace_hint(range.first, std::move(val));
        else
            belief_locations_[k].emplace(std::move(val));
    }
}

std::pair<eris_id_t, BookState> FileStorage::readBook() const {
    auto pair = std::make_pair<eris_id_t, BookState>(
            read_u64(), BookState(dimensions_));

    BookState &b = pair.second;
    b.id = pair.first;

    b.author = read_u64();
    for (uint32_t d = 0; d < dimensions_; d++) {
        b.position[d] = read_dbl();
    }

    b.quality = read_dbl();
    b.price = read_dbl();
    b.market = not std::isnan(b.price);
    b.revenue = read_dbl();
    b.revenue_lifetime = read_dbl();
    b.sales = read_u64();
    b.sales_lifetime = read_u64();
    b.copies = read_u64();
    b.age = read_u64();
    b.lifetime = read_u64();

    return pair;
}

void FileStorage::writeBook(const BookState &b) {
    write_u64(b.id);
    write_u64(b.author);

    for (size_t i = 0; i < b.position.dimensions; i++)
        write_dbl(b.position[i]);

    write_dbl(b.quality);
    write_dbl(b.price);
    write_dbl(b.revenue);
    write_dbl(b.revenue_lifetime);
    write_u64(b.sales);
    write_u64(b.sales_lifetime);
    write_u64(b.copies);
    write_u64(b.age);
    write_u64(b.lifetime);
}

}}
