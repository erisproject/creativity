#include "creativity/state/FileStorage.hpp"
#include "creativity/Creativity.hpp"
#include <cmath>
#include <sys/stat.h>

#ifdef ERIS_DEBUG
#define FILESTORAGE_DEBUG_WRITE_START \
    int64_t debug_start_location = f_.tellp();
#define FILESTORAGE_DEBUG_WRITE_CHECK(EXPECT) \
    int64_t debug_length = f_.tellp() - debug_start_location; \
    int64_t debug_expect = EXPECT; \
    if (debug_length != debug_expect) \
        throw std::runtime_error("Internal FileStorage error: " + std::string(__func__) + " wrote " + std::to_string(debug_length) + ", expected " + std::to_string(debug_expect));
#else
#define FILESTORAGE_DEBUG_WRITE_START
#define FILESTORAGE_DEBUG_WRITE_CHECK(EXPECT)
#endif

namespace creativity { namespace state {

using namespace eris;
using namespace creativity::belief;
using namespace Eigen;

/// Member definitions for header/cblock constants
constexpr unsigned int
        FileStorage::BLOCK_SIZE,
        FileStorage::HEADER::size,
        FileStorage::HEADER::states,
        FileStorage::CBLOCK::states,
        FileStorage::CBLOCK::next_cblock,
        FileStorage::LIBRARY::record_size,
        FileStorage::LIBRARY::block_records;
constexpr char
        FileStorage::HEADER::fileid[],
        FileStorage::HEADER::test_value[],
        FileStorage::ZERO_BLOCK[];
constexpr uint64_t FileStorage::HEADER::u64_test;
constexpr int64_t
        FileStorage::HEADER::pos::fileid,
        FileStorage::HEADER::pos::filever,
        FileStorage::HEADER::pos::test_value,
        FileStorage::HEADER::pos::num_states,
        FileStorage::HEADER::pos::dimensions,
        FileStorage::HEADER::pos::readers,
        FileStorage::HEADER::pos::boundary,
        FileStorage::HEADER::pos::book_distance_sd,
        FileStorage::HEADER::pos::book_quality_sd,
        FileStorage::HEADER::pos::reader_step_sd,
        FileStorage::HEADER::pos::reader_creation_shape,
        FileStorage::HEADER::pos::reader_creation_scale_min,
        FileStorage::HEADER::pos::reader_creation_scale_max,
        FileStorage::HEADER::pos::cost_fixed,
        FileStorage::HEADER::pos::cost_unit,
        FileStorage::HEADER::pos::cost_piracy,
        FileStorage::HEADER::pos::income,
        FileStorage::HEADER::pos::piracy_begins,
        FileStorage::HEADER::pos::piracy_link_proportion,
        FileStorage::HEADER::pos::prior_weight,
        FileStorage::HEADER::pos::prior_weight_piracy,
        FileStorage::HEADER::pos::init_prob_write,
        FileStorage::HEADER::pos::init_q_min,
        FileStorage::HEADER::pos::init_q_max,
        FileStorage::HEADER::pos::init_p_min,
        FileStorage::HEADER::pos::init_p_max,
        FileStorage::HEADER::pos::init_prob_keep,
        FileStorage::HEADER::pos::init_keep_price,
        FileStorage::HEADER::pos::state_first,
        FileStorage::HEADER::pos::state_last,
        FileStorage::HEADER::pos::continuation,
        FileStorage::HEADER::i64_test;
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
    try {
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

        if (parse) {
            parseMetadata();
            have_settings = true;
        }
        else {
            writeEmptyHeader();
        }
    }
    catch (std::ios_base::failure &c) {
        throw std::ios_base::failure("Unable to open " + filename + ": " + strerror(errno));
    }
}

size_t FileStorage::size() const {
#ifndef CREATIVITY_DISABLE_THREADED_STORAGE
    std::unique_lock<std::mutex> lock(f_mutex_);
#endif
    return state_pos_.size();
}

void FileStorage::thread_insert(std::shared_ptr<const State> &&state) {
#ifndef CREATIVITY_DISABLE_THREADED_STORAGE
    std::unique_lock<std::mutex> lock(f_mutex_);
#endif
    f_.seekp(0, f_.end);
    auto location = f_.tellp();
    if (location < HEADER::size) {
        throw std::runtime_error("File size is invalid, unable to write new State");
    }
    else if (location == HEADER::size) {
        // This is the first insertion, so we need to write the library pointer block
        writeLibraryPointerBlock(state->readers);
    }

    // We're now at the end of the file, which is where we want to put it.

    // Before actually writing the state, we need to make sure the reader library blocks are
    // populated.
    for (auto &r : state->readers) {
        updateLibrary(r.second);
    }

    // Write out the record before storing its location: if writing fails, there may be some
    // unreferences garbage (a partially written state) at the end of the file, but the file will
    // still be consistent: the garbage won't be referenced anywhere and so will just be ignored
    // if the file is reloaded.
    location = newBlock();
    writeState(*state);
    addStateLocation(location);

    // If dimensions and boundary are 0, update them.
    if (settings_.dimensions == 0 and state->dimensions != 0) {
        f_.seekp(HEADER::pos::dimensions);
        write_u32(state->dimensions);
        settings_.dimensions = state->dimensions;
    }
    if (settings_.boundary == 0 and state->boundary != 0) {
        f_.seekp(HEADER::pos::boundary);
        write_value(state->boundary);
        settings_.boundary = state->boundary;
    }
}


std::shared_ptr<const State> FileStorage::load(eris_time_t t) const {
#ifndef CREATIVITY_DISABLE_THREADED_STORAGE
    std::unique_lock<std::mutex> lock(f_mutex_);
#endif

    if (t >= state_pos_.size()) return std::shared_ptr<const State>();

    auto state_pos = state_pos_[t];
    f_.seekg(state_pos);

    return readState();
}

void FileStorage::flush() {
    // Let StorageBackend worry about finishing the thread:
    StorageBackend::flush();
#ifndef CREATIVITY_DISABLE_THREADED_STORAGE
    std::unique_lock<std::mutex> lock(f_mutex_);
#endif
    // Ensure the filehandle output is flushed to disk
    f_.flush();
}

void FileStorage::writeEmptyHeader() {
    f_.seekp(0);
    f_.write(HEADER::fileid, sizeof HEADER::fileid); // 'CrSt' file signature
    write_u32(1); // file version
    f_.write(HEADER::test_value, sizeof HEADER::test_value); // All test values squished into one

    // Write out 0s for everything else:
    f_.write(ZERO_BLOCK, HEADER::size - HEADER::pos::num_states);
}

void FileStorage::readSettings(CreativitySettings &settings) const {
    settings = settings_;
}

void FileStorage::writeSettings(const CreativitySettings &settings) {
    if (have_settings) {
        if (settings_.dimensions != settings.dimensions)
            throw std::logic_error("Cannot overwrite settings with a difference number of dimensions");
        if (settings_.boundary != settings.boundary)
            throw std::logic_error("Cannot overwrite settings with a different boundary");
    }

#ifndef CREATIVITY_DISABLE_THREADED_STORAGE
    std::unique_lock<std::mutex> lock(f_mutex_);
#endif

    f_.seekp(HEADER::pos::num_states);
    write_u32(state_pos_.size()); // Number of states
    write_value(settings.dimensions);
    write_value(settings.readers);
    write_value(settings.boundary);
    write_value(settings.book_distance_sd);
    write_value(settings.book_quality_sd);
    write_value(settings.reader_step_sd);
    write_value(settings.reader_creation_shape);
    write_value(settings.reader_creation_scale_min);
    write_value(settings.reader_creation_scale_max);
    write_value(settings.cost_fixed);
    write_value(settings.cost_unit);
    write_value(settings.cost_piracy);
    write_value(settings.income);
    write_value(settings.piracy_begins);
    write_value(settings.piracy_link_proportion);
    write_value(settings.prior_weight);
    write_value(settings.prior_weight_piracy);
    write_value(settings.initial.prob_write);
    write_value(settings.initial.q_min);
    write_value(settings.initial.q_max);
    write_value(settings.initial.p_min);
    write_value(settings.initial.p_max);
    write_value(settings.initial.prob_keep);
    write_value(settings.initial.keep_price);
    write_value(settings.initial.belief_threshold);

    // 4 bytes padding:
    write_u32(0); // Unused padding value

    if (f_.tellp() != HEADER::pos::state_first) {
        // If this exception occurs, something in the above sequence is wrong.
        throw std::runtime_error("Header writing failed: header size != " + std::to_string(HEADER::size) + " bytes");
    }

    // Copy the settings (so that readSettings() will return the right thing)
    settings_ = settings;

    // The rest is state positions, which we aren't supposed to touch.
}

void FileStorage::parseLibraryPointerBlock() {
    f_.seekg(HEADER::size, f_.beg);
    std::vector<std::pair<uint32_t, int64_t>> blocks;
    for (uint32_t i = 0; i < settings_.readers; i++) {
        auto rid = read_u32();
        auto loc = read_i64();
        blocks.emplace_back(std::move(rid), std::move(loc));
    }

    for (auto &b : blocks) {
        auto &rid = b.first;
        f_.seekg(b.second, f_.beg);
        lib_data data(0); // Will come back to reset the "0"
        while (true) {
            if (data.records_remaining > 0) {
                uint32_t bid = read_u32();
                if (bid > 0) {
                    uint32_t acquired = read_u32();
                    double quality = read_dbl();
                    uint8_t status = read_u8();
                    BookCopy copy(
                            quality,
                            status == 0 ? BookCopy::Status::wrote :  status == 1 ? BookCopy::Status::purchased :  BookCopy::Status::pirated,
                            acquired);
                    auto ins = data.library.emplace(std::move(bid), std::move(copy));
                    data.library_acq_sorted.emplace(std::ref(*ins.first));
                    data.records_remaining--;
                }
                else {
                    // The book id is 0, which means there are no more library books, so done.
                    // The next pos is the u32 (=0) we just read
                    data.pos_next = -4 + f_.tellg();
                    break;
                }
            }
            else {
                // Nothing left in this block, so the next thing is a pointer to the next block
                int64_t next_block = read_i64();
                if (next_block > 0) {
                    f_.seekg(next_block, f_.beg);
                    data.records_remaining = LIBRARY::block_records;
                }
                else {
                    // The pointer is 0, which means there is no next block, so done.
                    // the next pos is the block (=0) we just read
                    data.pos_next = -8 + f_.tellg();
                    break;
                }
            }
        }

        reader_lib_.emplace((unsigned int) rid, std::move(data));
    }
}

int64_t FileStorage::newBlock() {
    f_.seekp(0, f_.end);
    int64_t location = f_.tellp();
    int padding = location % BLOCK_SIZE;
    if (padding > 0) {
        f_.write(ZERO_BLOCK, padding);
        location += padding;
    }
    return location;
}

void FileStorage::writeLibraryPointerBlock(const std::unordered_map<eris_id_t, ReaderState> &readers) {
    f_.seekp(0, f_.end);
    if (f_.tellp() != HEADER::size)
        throwParseError("writing library pointer block failed: file is not just the header");
    // We're going to write out the first library blocks immediately after this section, so figure
    // out where the first one is: it'll be at the end, padding out to the next block.
    int64_t first_lib_block = HEADER::size + 12 * readers.size();
    first_lib_block += first_lib_block % BLOCK_SIZE;
    int64_t next_lib_block = first_lib_block;
    for (auto &r : readers) {
        write_u32(r.first);
        write_i64(next_lib_block);
        reader_lib_.emplace((unsigned int) r.first, lib_data(next_lib_block));
        next_lib_block += BLOCK_SIZE;
    }
    // Call newBlock() to add the padding, and check to make sure it's in the right place
    if (newBlock() != first_lib_block)
        throw std::runtime_error("writeLibraryPointBlock failed: next block isn't where it should be!");

    // Now write out the (empty) library blocks we just promised were there
    for (unsigned int i = 0; i < readers.size(); i++) {
        f_.write(ZERO_BLOCK, BLOCK_SIZE);
    }
}

void FileStorage::updateLibrary(const ReaderState &r) {
    lib_data *libd;
    try { libd = &reader_lib_.at(r.id); }
    catch (const std::out_of_range &e) {
        return throwParseError("Write failed: reader not found in library (perhaps readers aren't constant over sim periods?)");
    }
    auto &lib = *libd;
    for (auto &l : r.library) {
        bool first = true;
        if (lib.library.count(l.first) == 0) { // Need to insert
            if (lib.records_remaining == 0) {
                // There is no more space in the current block, so create a new one:
                int64_t next_block = newBlock();
                f_.write(ZERO_BLOCK, BLOCK_SIZE);
                // Now seek back to the end of the old record and write the new block location:
                f_.seekp(lib.pos_next);
                write_i64(next_block);
                // Now seek to the first location in the new block, to write the value (below)
                f_.seekp(next_block);
                lib.records_remaining = LIBRARY::block_records;
                lib.pos_next = next_block;
                first = false;
            }
            else if (first) {
                // If this is the first insertion for this reader, the current file pointer could be
                // anywhere; seek to the next positions.  (After the first insertion, the file
                // pointer will be just after the previous insertion, which is where we want it)
                f_.seekp(lib.pos_next);
                first = false;
            }
            // Otherwise we're already there from the last write
            write_u32(l.first);
            write_u32(l.second.acquired);
            write_dbl(l.second.quality);
            write_u8(l.second.wrote() ? 0 : l.second.purchased() ? 1 : 2);
            lib.records_remaining--;
            lib.pos_next += LIBRARY::record_size;

            auto ins = lib.library.emplace(l.first, l.second);
            lib.library_acq_sorted.emplace(std::ref(*ins.first));
        }
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
    write_u64(location);

    f_.seekp(HEADER::pos::num_states);
    state_pos_.push_back(location);
    write_u32(curr_states + 1);
}

void FileStorage::createContinuationBlock() {
    int64_t location = newBlock();
    f_.write(ZERO_BLOCK, BLOCK_SIZE);
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

    // number of states:
    auto num_states = parse_value<uint32_t>(block[HEADER::pos::num_states]);

    parse_value(block[HEADER::pos::dimensions], settings_.dimensions);
    parse_value(block[HEADER::pos::readers], settings_.readers);
    parse_value(block[HEADER::pos::boundary], settings_.boundary);
    parse_value(block[HEADER::pos::book_distance_sd], settings_.book_distance_sd);
    parse_value(block[HEADER::pos::book_quality_sd], settings_.book_quality_sd);
    parse_value(block[HEADER::pos::reader_step_sd], settings_.reader_step_sd);
    parse_value(block[HEADER::pos::reader_creation_shape], settings_.reader_creation_shape);
    parse_value(block[HEADER::pos::reader_creation_scale_min], settings_.reader_creation_scale_min);
    parse_value(block[HEADER::pos::reader_creation_scale_max], settings_.reader_creation_scale_max);
    parse_value(block[HEADER::pos::cost_fixed], settings_.cost_fixed);
    parse_value(block[HEADER::pos::cost_unit], settings_.cost_unit);
    parse_value(block[HEADER::pos::cost_piracy], settings_.cost_piracy);
    parse_value(block[HEADER::pos::income], settings_.income);
    parse_value(block[HEADER::pos::piracy_begins], settings_.piracy_begins);
    parse_value(block[HEADER::pos::prior_weight], settings_.prior_weight);
    parse_value(block[HEADER::pos::prior_weight_piracy], settings_.prior_weight_piracy);
    parse_value(block[HEADER::pos::piracy_link_proportion], settings_.piracy_link_proportion);
    parse_value(block[HEADER::pos::init_prob_write], settings_.initial.prob_write);
    parse_value(block[HEADER::pos::init_q_min], settings_.initial.q_min);
    parse_value(block[HEADER::pos::init_q_max], settings_.initial.q_max);
    parse_value(block[HEADER::pos::init_p_min], settings_.initial.p_min);
    parse_value(block[HEADER::pos::init_p_max], settings_.initial.p_max);
    parse_value(block[HEADER::pos::init_prob_keep], settings_.initial.prob_keep);
    parse_value(block[HEADER::pos::init_keep_price], settings_.initial.keep_price);
    parse_value(block[HEADER::pos::init_belief_threshold], settings_.initial.belief_threshold);

    if (settings_.dimensions == 0) throwParseError("found invalid dimensions == 0");
    if (settings_.readers == 0) throwParseError("found invalid readers == 0");
    if (settings_.boundary <= 0) throwParseError("found invalid (non-positive) boundary");
    if (settings_.book_distance_sd < 0) throwParseError("found invalid (negative) book_distance_sd");
    if (settings_.book_quality_sd < 0) throwParseError("found invalid (negative) book_quality_sd");

    if (num_states > 0) {
        // A reader library block is written immediately after the header before writing the first
        // state, so parse it.
        parseLibraryPointerBlock();
    }

    state_pos_.reserve(num_states);

    size_t header_states = std::min(num_states, HEADER::states);
    parseStateLocations(block[HEADER::pos::state_first], header_states, file_size);

    uint32_t remaining = num_states - header_states;
    while (remaining > 0) {
        auto cont = parse_value<std::streampos>(block[HEADER::pos::continuation]);
        if (cont < HEADER::size or cont >= file_size)
            throwParseError("found invalid continuation location");
        cont_pos_.push_back(cont);

        char cblock[BLOCK_SIZE];
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
    State *st_ptr = new State();
    State &state = *st_ptr;
    std::shared_ptr<const State> shst(st_ptr);

    state.t = read_u32();
    state.dimensions = settings_.dimensions;
    state.boundary = settings_.boundary;

    auto num_readers = read_u32();
    state.readers.reserve(num_readers);
    for (uint32_t i = 0; i < num_readers; i++) {
        state.readers.insert(readReader(state.t));
    }

    auto num_books = read_u32();
    state.books.reserve(num_books);
    for (uint32_t i = 0; i < num_books; i++) {
        state.books.insert(readBook());
    }

    return shst;
}

std::pair<eris_id_t, ReaderState> FileStorage::readReader(eris_time_t t) const {
    auto pair = std::make_pair<eris_id_t, ReaderState>(
            read_u32(),
            ReaderState(settings_.dimensions));

    ReaderState &r = pair.second;
    r.id = pair.first;

    // Position
    for (uint32_t d = 0; d < settings_.dimensions; d++)
        r.position[d] = read_dbl();

    // Friends
    auto num_friends = read_u32();
    r.friends.reserve(num_friends);
    for (uint32_t i = 0; i < num_friends; i++)
        r.friends.insert(read_u32());

    // Library
    auto &lib = reader_lib_.at(r.id);
    for (auto &l : lib.library_acq_sorted) {
        auto &id = l.get().first;
        auto &bc = l.get().second;
        if (bc.acquired > t) break;
        r.library.emplace(id, bc);
        if (bc.wrote())
            r.wrote.insert(id);
        else if (bc.acquired == t) // else if because new_books isn't supposed to contain self-authored books
            r.new_books.insert(id);
    }
    r.updateLibraryCounts(t);

    // Utility
    r.u = read_dbl();
    r.u_lifetime = read_dbl();
    // Costs
    r.cost_fixed = read_dbl();
    r.cost_unit = read_dbl();
    r.cost_piracy = read_dbl();
    r.income = read_dbl();
    r.creation_shape = read_dbl();
    r.creation_scale = read_dbl();

    // Beliefs
    belief_data belief = readBelief();
    if (belief.K > 0) {
        r.profit = belief.noninformative
            ? Profit(settings_.dimensions, belief.K)
            : Profit(settings_.dimensions, belief.beta, belief.s2, belief.V, belief.n);
        r.profit.draw_rejection_success = belief.draw_success_cumulative;
        r.profit.draw_rejection_discards = belief.draw_discards_cumulative;
    }

    belief = readBelief();
    if (belief.K > 0) {
        r.profit_extrap = belief.noninformative
            ? Profit(settings_.dimensions, belief.K)
            : Profit(settings_.dimensions, belief.beta, belief.s2, belief.V, belief.n);
        r.profit_extrap.draw_rejection_success = belief.draw_success_cumulative;
        r.profit_extrap.draw_rejection_discards = belief.draw_discards_cumulative;
    }

    belief = readBelief();
    if (belief.K > 0) {
        r.demand = belief.noninformative
            ? Demand(settings_.dimensions, belief.K)
            : Demand(settings_.dimensions, belief.beta, belief.s2, belief.V, belief.n);
        r.demand.draw_rejection_success = belief.draw_success_cumulative;
        r.demand.draw_rejection_discards = belief.draw_discards_cumulative;
    }

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
        r.profit_stream.emplace((unsigned int) belief.beta.rows(), belief.noninformative
                ? ProfitStream(belief.K)
                : ProfitStream(belief.beta, belief.s2, belief.V, belief.n));
    }

    return pair;
}

FileStorage::belief_data FileStorage::readBelief() const {
    belief_data belief;
    auto k = read_i8();
    // Handle magic values:
    // — 0 is not allowed.
    // — Positive values are the number of belief parameters.  We accept up to 127, though
    //   writeBelief() only writes beliefs with up to 120 parameters
    // — -1 through -120 are non-informative beliefs with 1 to 120 parameters
    // — -128 is a default-constructed belief object (i.e. not usable)
    // — Anything else (-121 to -127) is reserved for future use.
    if (k == -128) {
        belief.K = 0;
        return belief;
    }
    else if (k >= -120 and k < 0) {
        belief.K = (uint32_t) -k;
        belief.noninformative = true;
        return belief;
    }
    else if (k <= 0) {
        throwParseError("Found invalid belief record with k = " + std::to_string(k));
    }
    // Otherwise k is the right value.

    belief.K = k;
    belief.noninformative = false;

    // Next up is a status field
    uint8_t status = read_u8();
    // The lowest-value bit indicates this is a restricted belief
    // information)
    bool restricted_model = status & 1;
    bool last_draw_was_gibbs = status & 2;
    belief.fullyinformative = status & 4;

    // The first K elements are beta values
    belief.beta = VectorXd(k);
    for (unsigned int i = 0; i < belief.K; i++)
        belief.beta[i] = read_dbl();

    // Then s2 and n:
    belief.s2 = read_dbl();
    belief.n = read_dbl();

    // Then K*(K+1)/2 V values (but we set them symmetrically in V)
    belief.V = MatrixXd(k, k);
    for (unsigned int r = 0; r < belief.K; r++) {
        for (unsigned int c = 0; c <= r; c++) {
            double cov = read_dbl();
            belief.V(r,c) = cov;
            if (c != r) belief.V(c,r) = cov;
        }
    }

    // If the model was not fully informative, we immediately get a value r < K telling us how many
    // rows worth of independent data follow (i.e. r*K values follow).
    if (not belief.fullyinformative) {
        uint8_t indep_rows = read_u8();
        belief.indep_data = Matrix<double, Dynamic, Dynamic, RowMajor>(indep_rows, k);
        for (uint8_t i = 0; i < indep_rows; i++) for (unsigned int j = 0; j < belief.K; j++) {
            belief.indep_data(i, j) = read_dbl();
        }
    }

    // If this was a restricted model, we read the draw discards and success values
    if (restricted_model) {
        belief.draw_success_cumulative = read_u32();
        belief.draw_discards_cumulative = read_u32();
        belief.draw_gibbs = last_draw_was_gibbs;
    }

    return belief;
}

void FileStorage::writeState(const State &state) {

    // Now write the state:
    write_u32(state.t);
    if (state.readers.size() != settings_.readers)
        throw std::runtime_error("FileStorage error: cannot write a simulation where the number of readers changes");

    write_u32(state.readers.size());
    for (auto &r : state.readers) {
        writeReader(r.second);
    }
    write_u32(state.books.size());
    for (auto &b : state.books) {
        writeBook(b.second);
    }
}

void FileStorage::writeReader(const ReaderState &r) {
    FILESTORAGE_DEBUG_WRITE_START
    if (r.id > std::numeric_limits<uint32_t>::max())
        throw std::runtime_error("FileStorage error: cannot handle reader ids > 32 bits");
    write_u32(r.id);

    // Position
    for (size_t i = 0; i < r.position.dimensions; i++)
        write_dbl(r.position[i]);

    // Friends
    write_u32(r.friends.size());
    for (auto &f : r.friends)
        write_u32(f);

    // Utility
    write_value(r.u);
    write_value(r.u_lifetime);
    // Costs
    write_value(r.cost_fixed);
    write_value(r.cost_unit);
    write_value(r.cost_piracy);
    write_value(r.income);
    // Creation parameters
    write_value(r.creation_shape);
    write_value(r.creation_scale);

    FILESTORAGE_DEBUG_WRITE_CHECK(4*(1+2*settings_.dimensions+1+r.friends.size()+2+2+2+2+2+2+2+2))

    // Beliefs
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
    FILESTORAGE_DEBUG_WRITE_START
    auto &k = m.K();
    if (k > 120) {
        throw std::runtime_error("creativity::state::FileStorage cannot handle beliefs with K > 120");
    }
    else if (k == 0) {
        // If K = 0 (a default-constructed, non-model object), we just write out the special value and we're done
        write_i8(-128);
        return;
    }
    else if (m.noninformative()) {
        // Non-informative models are easy, too: just write out -K
        write_i8(-(int8_t) k);
        return;
    }

    // Otherwise write out the full record, starting with the positive size:
    write_i8(k);

    bool restricted_model = false;
    // First up is the status field.  Currently we have one status bit for a restricted model, and
    // one for the last draw mode
    const LinearRestricted *lr = dynamic_cast<const LinearRestricted*>(&m);
    if (lr) restricted_model = true;

    uint8_t status = 0;
    if (restricted_model) {
        status |= 1;
        if (lr->last_draw_mode == LinearRestricted::DrawMode::Gibbs) status |= 2;
    }
    if (m.fullyInformative()) {
        status |= 4;
    }

    // Status field
    write_u8(status);

    auto &beta = m.beta();
    // First K elements are the beta values
    for (unsigned int i = 0; i < k; i++)
        write_dbl(beta[i]);

    // Then s2 and n:
    write_dbl(m.s2());
    write_dbl(m.n());

    auto &V = m.V();
    // The last k*(k+1)/2 are the lower triangle of the V matrix, in row major order
    for (unsigned int r = 0; r < k; r++) {
        for (unsigned int c = 0; c <= r; c++) {
            write_dbl(V(r,c));
        }
    }

    if (not m.fullyInformative()) {
        auto indep_data = m.indepDataRows();
        write_u8(indep_data.rows());
        for (int i = 0; i < indep_data.rows(); i++) for (unsigned j = 0; j < k; j++) {
            write_dbl(indep_data(i,j));
        }
    }

    if (restricted_model) {
        write_u32(lr->draw_rejection_success);
        write_u32(lr->draw_rejection_discards);
    }

    FILESTORAGE_DEBUG_WRITE_CHECK(2+8*(2+k+k*(k+1)/2) + (restricted_model ? 4*2 : 0) +
            (m.fullyInformative() ? 0 : 1 + 8*m.indepDataRows().rows()*k))
}

std::pair<eris_id_t, BookState> FileStorage::readBook() const {
    auto pair = std::make_pair<eris_id_t, BookState>(
            read_u32(), BookState(settings_.dimensions));

    BookState &b = pair.second;
    b.id = pair.first;

    b.author = read_u32();
    for (uint32_t d = 0; d < settings_.dimensions; d++) {
        b.position[d] = read_dbl();
    }

    b.quality = read_dbl();
    b.price = read_dbl();
    b.revenue = read_dbl();
    b.revenue_lifetime = read_dbl();
    b.sales = read_u32();
    b.sales_lifetime = read_u32();
    b.pirated = read_u32();
    b.pirated_lifetime = read_u32();
    b.created = read_u32();
    b.lifetime = read_u32();

    return pair;
}

void FileStorage::writeBook(const BookState &b) {
    if (b.id > std::numeric_limits<uint32_t>::max())
        throw std::runtime_error("FileStorage error: cannot handle reader ids > 32 bits");
    write_u32(b.id);
    write_u32(b.author);

    for (size_t i = 0; i < b.position.dimensions; i++)
        write_dbl(b.position[i]);

    write_dbl(b.quality);
    write_dbl(b.price);
    write_dbl(b.revenue);
    write_dbl(b.revenue_lifetime);
    write_u32(b.sales);
    write_u32(b.sales_lifetime);
    write_u32(b.pirated);
    write_u32(b.pirated_lifetime);
    write_u32(b.created);
    write_u32(b.lifetime);
}

}}
