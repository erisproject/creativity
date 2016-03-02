#include "creativity/state/FileStorage.hpp"
#include <cerrno>
#include <cstring>
#include <limits>
#include <map>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <eris/Position.hpp>
#include "creativity/belief/Demand.hpp"
#include "creativity/belief/Profit.hpp"
#include "creativity/belief/ProfitStream.hpp"
#include <boost/math/constants/constants.hpp>
#include <boost/filesystem/operations.hpp>
extern "C" {
#include "lzma.h"
}

#ifdef ERIS_DEBUG
#define FILESTORAGE_DEBUG_WRITE_START \
    int64_t debug_start_location = f_->tellp();
#define FILESTORAGE_DEBUG_WRITE_CHECK(EXPECT) \
    int64_t debug_length = f_->tellp() - debug_start_location; \
    int64_t debug_expect = EXPECT; \
    if (debug_length != debug_expect) \
        throw std::runtime_error("Internal FileStorage error: " + std::string(__func__) + " wrote " + std::to_string(debug_length) + ", expected " + std::to_string(debug_expect));
#else
#define FILESTORAGE_DEBUG_WRITE_START
#define FILESTORAGE_DEBUG_WRITE_CHECK(EXPECT)
#endif

namespace creativity { namespace state {

using namespace eris;
using namespace eris::belief;
using namespace creativity::belief;
using namespace Eigen;
using DrawMode = BayesianLinearRestricted::DrawMode;
namespace fs = boost::filesystem;

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
constexpr  int64_t FileStorage::HEADER::i64_test;
constexpr uint32_t FileStorage::HEADER::u32_test;
constexpr  int32_t FileStorage::HEADER::i32_test;
constexpr double FileStorage::HEADER::dbl_test;
constexpr int64_t
        FileStorage::HEADER::pos::fileid,
        FileStorage::HEADER::pos::filever,
        FileStorage::HEADER::pos::test_value,
        FileStorage::HEADER::pos::num_states,
        FileStorage::HEADER::pos::dimensions,
        FileStorage::HEADER::pos::readers,
        FileStorage::HEADER::pos::boundary,
        FileStorage::HEADER::pos::book_distance_mean,
        FileStorage::HEADER::pos::book_quality_sd,
        FileStorage::HEADER::pos::reader_step_mean,
        FileStorage::HEADER::pos::reader_creation_shape,
        FileStorage::HEADER::pos::reader_creation_scale_min,
        FileStorage::HEADER::pos::reader_creation_scale_range,
        FileStorage::HEADER::pos::creation_time,
        FileStorage::HEADER::pos::cost_market,
        FileStorage::HEADER::pos::cost_unit,
        FileStorage::HEADER::pos::cost_piracy,
        FileStorage::HEADER::pos::income,
        FileStorage::HEADER::pos::piracy_begins,
        FileStorage::HEADER::pos::piracy_link_proportion,
        FileStorage::HEADER::pos::prior_scale,
        FileStorage::HEADER::pos::prior_scale_piracy,
        FileStorage::HEADER::pos::prior_scale_burnin,
        FileStorage::HEADER::pos::burnin_periods,
        FileStorage::HEADER::pos::init_prob_write,
        FileStorage::HEADER::pos::init_l_min,
        FileStorage::HEADER::pos::init_l_range,
        FileStorage::HEADER::pos::init_p_min,
        FileStorage::HEADER::pos::init_p_range,
        FileStorage::HEADER::pos::init_prob_keep,
        FileStorage::HEADER::pos::init_keep_price,
        FileStorage::HEADER::pos::public_sharing_begins,
        FileStorage::HEADER::pos::public_sharing_tax,
        FileStorage::HEADER::pos::prior_scale_public_sharing,
        FileStorage::HEADER::pos::creation_fixed,
        FileStorage::HEADER::pos::state_first,
        FileStorage::HEADER::pos::state_last,
        FileStorage::HEADER::pos::continuation;

constexpr std::ios_base::openmode FileStorage::open_readonly, FileStorage::open_overwrite, FileStorage::open_readwrite;

FileStorage::ParseError::ParseError(const std::string &message) : std::runtime_error(message) {}

void FileStorage::throwParseError(const std::string& message) const {
    decltype(f_->tellg()) pos;
    try { pos = f_->tellg(); }
    catch (std::ios_base::failure &f) { pos = -1; } // -1 is also be returned by .tellg() for pre-existing failures
    throw ParseError("Parsing file failed [pos=" + (pos == -1 ? std::string{"(error)"} : std::to_string(pos)) + "]: " + message);
}

FileStorage::FileStorage() {
    f_.reset(new std::stringstream(open_readwrite));
    initialize(true);
}

FileStorage::FileStorage(const std::string &filename, MODE mode) {
    try {
        auto *f = new std::fstream;
        f_.reset(f);
        f_->exceptions(f_->failbit | f_->badbit);

        bool empty = true;
        switch (mode) {
            case MODE::READONLY:
                f->open(filename, open_readonly);
                empty = false;
                break;
            case MODE::READ:
                f->open(filename, open_readwrite);
                empty = false;
                break;
            case MODE::OVERWRITE:
                f->open(filename, open_overwrite);
                empty = true;
                break;
            case MODE::APPEND:
                if (fs::exists(filename)) {
                    f->open(filename, open_readwrite);
                    f->seekg(0, f_->end);
                    empty = f->tellp() == 0;
                }
                else {
                    // Open in overwrite mode if the file doesn't exist (open in read-write mode
                    // fails if the file doesn't exist (unless also truncating)).
                    f->open(filename, open_overwrite);
                    empty = true;
                }
                break;
        }

        initialize(empty);
    }
    catch (std::ios_base::failure &c) {
        throw std::ios_base::failure("Unable to open " + filename + ": " + strerror(errno));
    }
}

FileStorage::FileStorage(const std::string &filename, bool xz_to_ram, bool cr_to_ram, const std::string &tmpfile) {
    std::string erroneous_file = filename;
    try {
        std::unique_ptr<std::fstream> f(new std::fstream);
        f->exceptions(f->failbit | f->badbit);
        f->open(filename, open_readonly);
        char magic[6];
        // Make sure we fail to compile if fileid ever changes to something longer than 6 bytes,
        // this code will break:
        static_assert(sizeof(HEADER::fileid) == 4, "Internal error: expected .crstate header size != 4");

        // A .crstate starts with CrSt; an xz-compressed file starts with the 6 bytes: 0xFD "7zXZ" 0x0
        // If we eof here, we'll throw: good.
        f->read(magic, 6);
        f->seekg(0);
        bool is_crst = true;
        for (size_t i = 0; i < sizeof(HEADER::fileid); i++) {
            if (magic[i] != HEADER::fileid[i]) {
                is_crst = false;
                break;
            }
        }
        if (not is_crst and not (magic[0] == (char) 0xFD and magic[1] == '7' and magic[2] == 'z' and magic[3] == 'X' and magic[4] == 'Z' and magic[5] == 0)) {
            throw std::runtime_error("Unable to parse " + filename + ": neither .crstate nor .xz file signature found");
        }

        if (is_crst) {
            if (cr_to_ram) {
                f_.reset(new std::stringstream(open_readwrite));
                f_->exceptions(f_->failbit | f_->badbit);
                *f_ << f->rdbuf();
            }
            else {
                f_ = std::move(f);
            }
        }
        else {
            // XZ compressed
            if (xz_to_ram) {
                f_.reset(new std::stringstream(open_readwrite));
                f_->exceptions(f_->failbit | f_->badbit);
            }
            else {
                auto tmpf = new std::fstream;
                tmpf->exceptions(f_->failbit | f_->badbit);
                tmpfile_path_ = tmpfile.empty()
                    ? (fs::temp_directory_path() / fs::unique_path("creativity-%%%%-%%%%-%%%%-%%%%.crstate")).native()
                    : tmpfile;
                erroneous_file = tmpfile_path_;
                tmpf->open(tmpfile_path_, open_overwrite);
                f_.reset(tmpf);
            }
            decompressXZ(*f, *f_);
        }

        initialize(false);
    }
    catch (std::ios_base::failure &c) {
        throw std::ios_base::failure("Error opening " + erroneous_file + ": " + strerror(errno));
    }
}

FileStorage::~FileStorage() {
    if (f_ and not tmpfile_path_.empty()) {
        f_.reset();
        std::remove(tmpfile_path_.c_str());
    }
}

void FileStorage::initialize(bool empty) {
    if (empty)
        writeEmptyHeader();
    else
        parseMetadata();
}

size_t FileStorage::size() const {
    std::unique_lock<std::mutex> lock(f_mutex_);
    return state_pos_.size();
}

void FileStorage::thread_insert(std::shared_ptr<const State> &&state) {
    std::unique_lock<std::mutex> lock(f_mutex_);
    f_->seekp(0, f_->end);
    auto location = f_->tellp();
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
        f_->seekp(HEADER::pos::dimensions);
        write_u32(state->dimensions);
        settings_.dimensions = state->dimensions;
    }
    if (settings_.boundary == 0 and state->boundary != 0) {
        f_->seekp(HEADER::pos::boundary);
        write_value(state->boundary);
        settings_.boundary = state->boundary;
    }
}


std::shared_ptr<const State> FileStorage::load(eris_time_t t) const {
    std::unique_lock<std::mutex> lock(f_mutex_);

    if (t >= state_pos_.size()) return std::shared_ptr<const State>();

    auto state_pos = state_pos_[t];
    f_->seekg(state_pos);

    return readState();
}

void FileStorage::storage_flush() {
    std::unique_lock<std::mutex> lock(f_mutex_);
    // Ensure the filehandle output is flushed to disk
    f_->flush();
}

void FileStorage::writeEmptyHeader() {
    f_->seekp(0);
    f_->write(HEADER::fileid, sizeof HEADER::fileid); // 'CrSt' file signature
    write_u32(version_); // file version
    f_->write(HEADER::test_value, sizeof HEADER::test_value); // All test values squished into one

    // Write out 0s for everything else:
    f_->write(ZERO_BLOCK, HEADER::size - HEADER::pos::num_states);
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

    std::unique_lock<std::mutex> lock(f_mutex_);

    f_->seekp(HEADER::pos::num_states);
    write_u32(state_pos_.size()); // Number of states
    write_value(settings.dimensions);
    write_value(settings.readers);
    write_value(settings.boundary);
    write_value(settings.book_distance_mean);
    write_value(settings.book_quality_sd);
    write_value(settings.reader_step_mean);
    write_value(settings.reader_creation_shape);
    write_value(settings.reader_creation_scale_min);
    write_value(settings.reader_creation_scale_range);
    write_value(settings.creation_time);
    write_value(settings.cost_market);
    write_value(settings.cost_unit);
    write_value(settings.cost_piracy);
    write_value(settings.income);
    write_value(settings.piracy_begins);
    write_value(settings.piracy_link_proportion);
    write_value(settings.prior_scale);
    write_value(settings.prior_scale_piracy);
    write_value(settings.prior_scale_burnin);
    write_value(settings.burnin_periods);
    write_value(settings.initial.prob_write);
    write_value(settings.initial.l_min);
    write_value(settings.initial.l_range);
    write_value(settings.initial.p_min);
    write_value(settings.initial.p_range);
    write_value(settings.initial.prob_keep);
    write_value(settings.initial.keep_price);
    write_value(settings.initial.belief_threshold);
    write_value(settings.public_sharing_begins);
    write_value(settings.public_sharing_tax);
    write_value(settings.prior_scale_public_sharing);
    write_value(settings.creation_fixed);

    // Uncommented when padding needed:
    //write_u32(0); // Unused padding value

    int64_t expect = HEADER::pos::state_first;
    if (f_->tellp() != expect) {
        // If this exception occurs, something in the above sequence is wrong.
        throw std::runtime_error("Header writing failed: header parameter block != " + std::to_string(expect) + " bytes");
    }

    // Copy the settings (so that readSettings() will return the right thing)
    settings_ = settings;

    // The rest is state positions, which we aren't supposed to touch.
}

void FileStorage::parseLibraryPointerBlock() {
    f_->seekg(HEADER::size, f_->beg);
    std::vector<std::pair<uint32_t, int64_t>> blocks;
    for (uint32_t i = 0; i < settings_.readers; i++) {
        auto rid = read_u32();
        auto loc = read_i64();
        blocks.emplace_back(std::move(rid), std::move(loc));
    }

    for (auto &b : blocks) {
        auto &rid = b.first;
        f_->seekg(b.second, f_->beg);
        lib_data data(0); // Will come back to reset the "0"
        while (true) {
            if (data.records_remaining > 0) {
                uint32_t bid = read_u32();
                if (bid > 0) {
                    uint32_t acquired = read_u32();
                    double quality = read_dbl();
                    uint8_t status = read_u8();
                    if (status > 3) throw std::runtime_error("Found invalid library status value `" + std::to_string(status) + "'");
                    BookCopy copy(
                            quality,
                            (   status == 0 ? BookCopy::Status::wrote :
                                status == 1 ? BookCopy::Status::purchased_market :
                                status == 2 ? BookCopy::Status::pirated :
                                BookCopy::Status::purchased_public),
                            acquired);
                    auto ins = data.library.emplace(std::move(bid), std::move(copy));
                    data.library_acq_sorted.emplace(std::ref(*ins.first));
                    data.records_remaining--;
                }
                else {
                    // The book id is 0, which means there are no more library books, so done.
                    // The next pos is the u32 (=0) we just read
                    data.pos_next = -4 + f_->tellg();
                    break;
                }
            }
            else {
                // Nothing left in this block, so the next thing is a pointer to the next block
                int64_t next_block = read_i64();
                if (next_block > 0) {
                    f_->seekg(next_block, f_->beg);
                    data.records_remaining = LIBRARY::block_records;
                }
                else {
                    // The pointer is 0, which means there is no next block, so done.
                    // the next pos is the block (=0) we just read
                    data.pos_next = -8 + f_->tellg();
                    break;
                }
            }
        }

        reader_lib_.emplace((unsigned int) rid, std::move(data));
    }
}

int64_t FileStorage::newBlock() {
    f_->seekp(0, f_->end);
    int64_t location = f_->tellp();
    int padding = location % BLOCK_SIZE;
    if (padding > 0) {
        f_->write(ZERO_BLOCK, padding);
        location += padding;
    }
    return location;
}

void FileStorage::writeLibraryPointerBlock(const std::map<eris_id_t, ReaderState> &readers) {
    f_->seekp(0, f_->end);
    if (f_->tellp() != HEADER::size)
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
        f_->write(ZERO_BLOCK, BLOCK_SIZE);
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
                f_->write(ZERO_BLOCK, BLOCK_SIZE);
                // Now seek back to the end of the old record and write the new block location:
                f_->seekp(lib.pos_next);
                write_i64(next_block);
                // Now seek to the first location in the new block, to write the value (below)
                f_->seekp(next_block);
                lib.records_remaining = LIBRARY::block_records;
                lib.pos_next = next_block;
                first = false;
            }
            else if (first) {
                // If this is the first insertion for this reader, the current file pointer could be
                // anywhere; seek to the next positions.  (After the first insertion, the file
                // pointer will be just after the previous insertion, which is where we want it)
                f_->seekp(lib.pos_next);
                first = false;
            }
            // Otherwise we're already there from the last write
            write_u32(l.first);
            write_u32(l.second.acquired);
            write_dbl(l.second.quality);
            write_u8(l.second.wrote() ? 0 : l.second.purchased_market() ? 1 : l.second.pirated() ? 2 : /* l.second.purchased_public() */ 3);
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
        f_->seekp(HEADER::pos::state_first + 8 * (int64_t) curr_states);
    }
    else {
        f_->seekp(cont_pos_.back() + (std::streampos) 8 * (past_header % CBLOCK::states));
    }
    write_u64(location);

    f_->seekp(HEADER::pos::num_states);
    state_pos_.push_back(location);
    write_u32(curr_states + 1);
}

void FileStorage::createContinuationBlock() {
    int64_t location = newBlock();
    f_->write(ZERO_BLOCK, BLOCK_SIZE);
    // Write the new block location either in the header (if this is the first cblock) or in the
    // previous cblock
    f_->seekp(cont_pos_.empty() ? HEADER::pos::continuation : (int64_t) cont_pos_.back() + CBLOCK::next_cblock);
    write_i64(location);
    cont_pos_.push_back(location);
}

void FileStorage::parseMetadata() {
    f_->seekg(0, f_->end);
    auto file_size = f_->tellg();
    if (file_size == 0)
        throwParseError("file contains no data");
    else if (file_size < HEADER::size)
        throwParseError("file is too small to contain header data");

    char block[HEADER::size];
    f_->seekg(0);
    f_->read(block, HEADER::size);

    // Make sure the CrSt file signature is found:
    for (unsigned int i = 0; i < 4; i++) {
        if (block[HEADER::pos::fileid+i] != HEADER::fileid[i])
            throwParseError("'CrSt' file signature not found");
    }

    auto version = parse_value<uint32_t>(block[HEADER::pos::filever]);
    if (version != version_)
        throwParseError("encountered unknown/unsupported file version " + std::to_string(version) + " (expected " + std::to_string(version_) + ")");

    // Okay, so we've got a CrSt version 2 file.  The version number lets the file format change
    // later, if necessary, in which case this code will need to handle the different versions.  v1
    // to v2 introduced incompatible changes, and so only v2 versions are handled here.  Future file
    // versions could conceivably handle multiple versions at once, provided the features are
    // comparable.

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

#define PARSE_VALUE(FIELD) parse_value(block[HEADER::pos::FIELD], settings_.FIELD)
#define PARSE_VALUE_INIT(FIELD) parse_value(block[HEADER::pos::init_##FIELD], settings_.initial.FIELD)
    PARSE_VALUE(dimensions);
    PARSE_VALUE(readers);
    PARSE_VALUE(boundary);
    PARSE_VALUE(book_distance_mean);
    PARSE_VALUE(book_quality_sd);
    PARSE_VALUE(reader_step_mean);
    PARSE_VALUE(reader_creation_shape);
    PARSE_VALUE(reader_creation_scale_min);
    PARSE_VALUE(reader_creation_scale_range);
    PARSE_VALUE(creation_time);
    PARSE_VALUE(cost_market);
    PARSE_VALUE(cost_unit);
    PARSE_VALUE(cost_piracy);
    PARSE_VALUE(income);
    PARSE_VALUE(piracy_begins);
    PARSE_VALUE(piracy_link_proportion);
    PARSE_VALUE(prior_scale);
    PARSE_VALUE(prior_scale_piracy);
    PARSE_VALUE(prior_scale_burnin);
    PARSE_VALUE(burnin_periods);
    PARSE_VALUE_INIT(prob_write);
    PARSE_VALUE_INIT(l_min);
    PARSE_VALUE_INIT(l_range);
    PARSE_VALUE_INIT(p_min);
    PARSE_VALUE_INIT(p_range);
    PARSE_VALUE_INIT(prob_keep);
    PARSE_VALUE_INIT(keep_price);
    PARSE_VALUE_INIT(belief_threshold);
    PARSE_VALUE(public_sharing_begins);
    PARSE_VALUE(public_sharing_tax);
    PARSE_VALUE(prior_scale_public_sharing);
    PARSE_VALUE(creation_fixed);
#undef PARSE_VALUE
#undef PARSE_VALUE_INIT

    if (settings_.dimensions == 0) throwParseError("found invalid dimensions == 0");
    if (settings_.readers == 0) throwParseError("found invalid readers == 0");
    if (settings_.boundary <= 0) throwParseError("found invalid (non-positive) boundary");
    if (settings_.book_distance_mean < 0) throwParseError("found invalid (negative) book_distance_mean");
    if (settings_.book_quality_sd < 0) throwParseError("found invalid (negative) book_quality_sd");
    if (settings_.reader_step_mean < 0) throwParseError("found invalid (negative) reader_step_mean");

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
        f_->seekg(cont, f_->beg);
        f_->read(cblock, sizeof cblock);

        uint32_t block_locations = std::min(remaining, CBLOCK::states);
        parseStateLocations(cblock[0], block_locations, file_size);

        remaining -= block_locations;
    }

    if (num_states != state_pos_.size())
        throwParseError("found " + std::to_string(state_pos_.size()) +
                " state data locations but expected " + std::to_string(num_states));

    have_settings = true;
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

    for (uint8_t type = read_u8(); type != TYPE_DONE; type = read_u8()) {
        switch (type) {
            case TYPE_READERS:
                {
                    auto num = read_u32();
                    for (uint32_t i = 0; i < num; i++) state.readers.insert(readReader(state.t));
                    break;
                }
            case TYPE_BOOKS:
                {
                    auto num = read_u32();
                    for (uint32_t i = 0; i < num; i++) state.books.insert(readBook());
                    break;
                }
            case TYPE_PUBLIC_TRACKER:
                {
                    if (state.publictracker) throw std::runtime_error("File error: found multiple public trackers");
                    state.publictracker = readPublicTracker();
                    break;
                }
            default:
                throw std::runtime_error("File error: found invalid/unsupported record type `" + std::to_string(type) + "'");
        }
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
    for (unsigned d = 0; d < settings_.dimensions; d++)
        r.position[d] = read_dbl();

    // Friends
    auto num_friends = read_u32();
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
    r.creation_shape = read_dbl();
    r.creation_scale = read_dbl();

    // Beliefs
    belief_data belief = readBelief();
    if (belief.K > 0) {
        r.profit = belief.noninformative
            ? Profit(belief.K)
            : Profit(belief.beta, belief.s2, belief.Vinv, belief.n);
        r.profit.last_draw_mode = belief.last_draw_mode;
        r.profit.draw_rejection_success = belief.draw_success_cumulative;
        r.profit.draw_rejection_discards = belief.draw_discards_cumulative;
    }

    belief = readBelief();
    if (belief.K > 0) {
        r.profit_extrap = belief.noninformative
            ? Profit(belief.K)
            : Profit(belief.beta, belief.s2, belief.Vinv, belief.n);
        r.profit_extrap.last_draw_mode = belief.last_draw_mode;
        r.profit_extrap.draw_rejection_success = belief.draw_success_cumulative;
        r.profit_extrap.draw_rejection_discards = belief.draw_discards_cumulative;
    }

    belief = readBelief();
    if (belief.K > 0) {
        r.demand = belief.noninformative
            ? Demand(belief.K)
            : Demand(belief.beta, belief.s2, belief.Vinv, belief.n);
        r.demand.last_draw_mode = belief.last_draw_mode;
        r.demand.draw_rejection_success = belief.draw_success_cumulative;
        r.demand.draw_rejection_discards = belief.draw_discards_cumulative;
    }

    auto pstream_locs = read_u32();
    for (uint32_t i = 0; i < pstream_locs; i++) {
        belief = readBelief();
        if (belief.K == 0) throwParseError("found illegal profitStream belief with K = 0 (i.e. default constructed model)");
        else if (r.profit_stream.count(belief.K)) throwParseError("found duplicate K value in profit_stream beliefs");
        r.profit_stream.emplace((unsigned int) belief.beta.rows(), belief.noninformative
                ? ProfitStream(belief.K)
                : ProfitStream(belief.beta, belief.s2, belief.Vinv, belief.n));
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
    else if (k <= 0) {
        throwParseError("Found invalid belief record with k = " + std::to_string(k));
    }
    // Otherwise k is the right value.

    belief.K = k;

    // Next up is a status field
    uint8_t status = read_u8();
    // The lowest-value bit indicates this is a restricted belief
    // information)
    bool restricted_model = status & 1;
    belief.noninformative = status & 4;
    belief.last_draw_mode =
        (restricted_model) ?
            (status & 2) ? DrawMode::Gibbs :
            (status & 8) ? DrawMode::Rejection :
            DrawMode::Auto
        : DrawMode::Auto;

    if (not belief.noninformative) {
        // The first K elements are beta values
        belief.beta = VectorXd(k);
        for (unsigned int i = 0; i < belief.K; i++)
            belief.beta[i] = read_dbl();

        // Then s2 and n:
        belief.s2 = read_dbl();
        belief.n = read_dbl();

        // Then K*(K+1)/2 V values (but we set them symmetrically in V)
        belief.Vinv = MatrixXd(k, k);
        for (unsigned int r = 0; r < belief.K; r++) {
            for (unsigned int c = 0; c <= r; c++) {
                double cov = read_dbl();
                belief.Vinv(r,c) = cov;
                if (c != r) belief.Vinv(c,r) = cov;
            }
        }

        // If this was a restricted model, we read the draw discards and success values
        if (restricted_model) {
            belief.draw_success_cumulative = read_u32();
            belief.draw_discards_cumulative = read_u32();
        }
    }

    return belief;
}

void FileStorage::writeState(const State &state) {
    if (state.readers.size() != settings_.readers)
        throw std::runtime_error("FileStorage error: cannot write a simulation where the number of readers changes");

    // First the state's `t`:
    write_u32(state.t);

    // NB: The order of these doesn't technically matter, but state.readers is an ordered map so
    // that we order by eris_id_t and thus get reproducible files for identical state data.
    if (not state.readers.empty()) {
        write_u8(TYPE_READERS);
        write_u32(state.readers.size());
        for (auto &r : state.readers) writeReader(r.second);
    }

    if (not state.books.empty()) {
        write_u8(TYPE_BOOKS);
        write_u32(state.books.size());
        for (auto &b : state.books) writeBook(b.second);
    }

    if (state.publictracker and state.publictracker->id != 0) {
        write_u8(TYPE_PUBLIC_TRACKER);
        writePublicTracker(*state.publictracker);
    }

    write_u8(TYPE_DONE);
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
    // Creation parameters
    write_value(r.creation_shape);
    write_value(r.creation_scale);

    FILESTORAGE_DEBUG_WRITE_CHECK(4+8*settings_.dimensions+4+4*r.friends.size()+8+8+8+8)

    // Beliefs
    writeBelief(r.profit);
    writeBelief(r.profit_extrap);
    writeBelief(r.demand);

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

void FileStorage::writeBelief(const BayesianLinear &m) {
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

    // Otherwise write out the full record, starting with the positive size:
    write_i8(k);

    bool restricted_model = false;
    // First up is the status field.  Currently we have:
    // 1 - restricted
    // 2 - (if restricted) last draw was gibbs
    // 4 - noninformative
    // 8 - (if restricted) last draw was rejection sampling
    const BayesianLinearRestricted *lr = dynamic_cast<const BayesianLinearRestricted*>(&m);
    if (lr) restricted_model = true;

    uint8_t status = 0;
    if (restricted_model) {
        status |= 1;
        if (lr->last_draw_mode == DrawMode::Gibbs) status |= 2;
        else if (lr->last_draw_mode == DrawMode::Rejection) status |= 8;
    }
    if (m.noninformative()) {
        status |= 4;
    }

    // Status field
    write_u8(status);

    if (not m.noninformative()) {
        auto &beta = m.beta();
        // First K elements are the beta values
        for (unsigned int i = 0; i < k; i++)
            write_dbl(beta[i]);

        // Then s2 and n:
        write_dbl(m.s2());
        write_dbl(m.n());

        const auto &Vinv = m.Vinv();
        // The last k*(k+1)/2 are the lower triangle of the V matrix, in row major order
        for (unsigned int r = 0; r < k; r++) {
            for (unsigned int c = 0; c <= r; c++) {
                write_dbl(Vinv(r,c));
            }
        }

        if (restricted_model) {
            write_u32(lr->draw_rejection_success);
            write_u32(lr->draw_rejection_discards);
        }
    }

    FILESTORAGE_DEBUG_WRITE_CHECK(
            2 + (
                m.noninformative()
                ? 0 : (8*(2+k+k*(k+1)/2) + (restricted_model ? 4*2 : 0))
                )
            );
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
    auto stat = read_u8();
    b.market_private = stat & 1;
    stat &= ~1;
    if (stat != 0) throw std::runtime_error("FileStorage error: found invalid/unknown book status bit");
    b.price = read_dbl();
    b.revenue = read_dbl();
    b.revenue_lifetime = read_dbl();
    b.prize = read_dbl();
    b.prize_lifetime = read_dbl();
    b.sales = read_u32();
    b.sales_lifetime_private = read_u32();
    b.sales_lifetime_public = read_u32();
    b.pirated = read_u32();
    b.pirated_lifetime = read_u32();
    b.created = read_u32();
    b.lifetime_private = read_u32();

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
    uint8_t status = 0;
    if (b.market_private) status |= 1;
    write_u8(status);
    write_dbl(b.price);
    write_dbl(b.revenue);
    write_dbl(b.revenue_lifetime);
    write_dbl(b.prize);
    write_dbl(b.prize_lifetime);
    write_u32(b.sales);
    write_u32(b.sales_lifetime_private);
    write_u32(b.sales_lifetime_public);
    write_u32(b.pirated);
    write_u32(b.pirated_lifetime);
    write_u32(b.created);
    write_u32(b.lifetime_private);
}

std::unique_ptr<PublicTrackerState> FileStorage::readPublicTracker() const {
    std::unique_ptr<PublicTrackerState> pts(new PublicTrackerState);
    pts->id = read_u32();
    pts->tax = read_dbl();
    pts->unspent = read_dbl();
    return pts;
}

void FileStorage::writePublicTracker(const PublicTrackerState &pt) {
    if (pt.id > std::numeric_limits<uint32_t>::max())
        throw std::runtime_error("FileStorage error: cannot handle eris ids > 32 bits");
    write_u32(pt.id);
    write_dbl(pt.tax);
    write_dbl(pt.unspent);
}

void FileStorage::copyTo(const std::string &filename) {
    std::unique_lock<std::mutex> lock(f_mutex_);
    std::ofstream out;
    out.exceptions(out.failbit | out.badbit);
    out.open(filename, open_overwrite);
    f_->seekg(0, f_->beg);
    out << f_->rdbuf();
}

void FileStorage::copyToXZ(const std::string &filename, uint32_t level) {
    std::unique_lock<std::mutex> lock(f_mutex_);
    std::ofstream xzout;
    xzout.exceptions(xzout.failbit | xzout.badbit);
    xzout.open(filename, open_overwrite);
    f_->seekg(0, f_->beg);
    compressXZ(*f_, xzout, level);
}

void FileStorage::compressXZ(std::istream &in, std::ostream &out, uint32_t level) {
    lzma_stream strm = LZMA_STREAM_INIT;
    // Some testing with .crstate files convinced me that -3 is optimal: it's quite fast (compared
    // to -4 and above), and the higher numbers offer only a couple extra percentage of compression
    // (and actually, -4 did worse).
    lzma_ret ret = lzma_easy_encoder(&strm, level, LZMA_CHECK_CRC64);

    if (ret != LZMA_OK)
        throw std::runtime_error(std::string("liblzma initialization failed: ") +
                (ret == LZMA_MEM_ERROR ? "Memory allocation failed" :
                 ret == LZMA_OPTIONS_ERROR ? "Specified compression level is invalid/not supported" :
                 ret == LZMA_UNSUPPORTED_CHECK ? "Specified integrity check is not supported" :
                 "An unknown error occured"));

    processXZ(&strm, &ret, in, out);
}

void FileStorage::processXZ(void *lzma_stream_ptr, void *lzma_ret_ptr, std::istream &in, std::ostream &out) {
    // The signature really take these references directly, but that makes liblzma-dev a
    // compile-time dependency of anything including FileStorage, while this way avoids that.
    // (Normally a forward declaration would take care of that, but I can't seem to forward-declare
    // C structs in C++ code).
    lzma_stream &strm = *static_cast<lzma_stream*>(lzma_stream_ptr);
    lzma_ret &ret = *static_cast<lzma_ret*>(lzma_ret_ptr);
    uint8_t inbuf[BUFSIZ], outbuf[BUFSIZ];
    strm.next_in = nullptr;
    strm.avail_in = 0;
    strm.next_out = outbuf;
    strm.avail_out = sizeof(outbuf);
    lzma_action action = LZMA_RUN;

    auto save_in_ex = in.exceptions(), save_out_ex = out.exceptions();
    // We don't want failbit in the input exceptions, because we want to be able to hit eof without
    // throwing an exception.
    in.exceptions(in.badbit);
    out.exceptions(out.badbit | out.failbit);

    try {
        while (true) {
            if (strm.avail_in == 0 and not in.eof()) {
                in.read(reinterpret_cast<char*>(inbuf), sizeof(inbuf));
                strm.next_in = inbuf;
                strm.avail_in = in.gcount();
                if (in.eof())
                    action = LZMA_FINISH;
            }

            ret = lzma_code(&strm, action);

            if (strm.avail_out == 0 or ret == LZMA_STREAM_END) {
                auto write_size = sizeof(outbuf) - strm.avail_out;
                out.write(reinterpret_cast<char*>(outbuf), write_size);

                strm.next_out = outbuf;
                strm.avail_out = sizeof(outbuf);
            }

            if (ret != LZMA_OK) {
                if (ret == LZMA_STREAM_END) break;

                throw std::runtime_error(std::string("liblzma compression/decompression failed: ") +
                        (ret == LZMA_MEM_ERROR ? "Memory allocation failed" :
                         ret == LZMA_DATA_ERROR ? "File size limits exceeded" :
			             ret == LZMA_FORMAT_ERROR ? "Input not in .xz format" :
			             ret == LZMA_OPTIONS_ERROR ? "Unsupported compression options" :
			             ret == LZMA_DATA_ERROR ? "Compressed file is corrupt" :
			             ret == LZMA_BUF_ERROR ? "Compressed file is truncated or otherwise corrupt" :
                         "An unknown error occured"));
            }
        }
    }
    catch (const std::ios_base::failure&) {
        if (not in.bad()) in.clear();
        in.exceptions(save_in_ex);
        out.exceptions(save_out_ex);
        throw;
    }
    in.clear();
    in.exceptions(save_in_ex);
    out.exceptions(save_out_ex);
}

void FileStorage::decompressXZ(std::istream &in, std::ostream &out) {
    lzma_stream strm = LZMA_STREAM_INIT;
    lzma_ret ret = lzma_stream_decoder(&strm, UINT64_MAX, LZMA_CONCATENATED);

    if (ret != LZMA_OK)
        throw std::runtime_error(std::string("liblzma initialization failed: ") +
                (ret == LZMA_MEM_ERROR ? "Memory allocation failed" :
                 ret == LZMA_OPTIONS_ERROR ? "Unsupported decompression flags" :
                 "An unknown error occured"));
    processXZ(&strm, &ret, in, out);
}

}}
