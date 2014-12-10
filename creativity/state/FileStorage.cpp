#include "creativity/state/FileStorage.hpp"
#include "creativity/Creativity.hpp"
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

FileStorage::FileStorage(const std::string &filename, MODE mode, CreativitySettings &settings_in) {
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
        settings_in = settings;
    }
    else {
        settings_ = settings_in;
        writeEmptyHeader();
    }
}

size_t FileStorage::size() const {
    return state_pos_.size();
}

void FileStorage::push_back(std::shared_ptr<const State> state) {
    // Make sure the state agrees with the states we already have
    if (settings.dimensions != 0 and state->dimensions != settings.dimensions)
        throw std::runtime_error("Cannot add state: state dimensions differ from storage dimensions");
    if (settings.boundary != 0 and state->boundary != settings.boundary)
        throw std::runtime_error("Cannot add state: state boundary differs from storage boundary");
    if (settings.readers != 0 and settings.readers != state->readers.size())
        throw std::runtime_error("Cannot add state: state #readers differs from storage #readers");

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

    if (settings.dimensions == 0 and state->dimensions != 0) {
        f_.seekp(HEADER::pos::dimensions);
        write_u32(state->dimensions);
        settings_.dimensions = state->dimensions;
    }
    if (settings.boundary == 0 and state->boundary != 0) {
        f_.seekp(HEADER::pos::boundary);
        write_value(state->boundary);
        settings_.boundary = state->boundary;
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

    write_u32(0); // Number of states
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

    // Currently no padding needed:
    //write_u32(0); // Unused padding value

    // State addresses and continuation location (all 0):
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

    // number of states:
    auto num_states = parse_value<uint32_t>(block[HEADER::pos::states]);

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

    if (settings.dimensions == 0) throwParseError("found invalid dimensions == 0");
    if (settings.readers == 0) throwParseError("found invalid readers == 0");
    if (settings.boundary <= 0) throwParseError("found invalid (non-positive) boundary");
    if (settings.book_distance_sd < 0) throwParseError("found invalid (negative) book_distance_sd");
    if (settings.book_quality_sd < 0) throwParseError("found invalid (negative) book_quality_sd");

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
    State *st_ptr = new State();
    State &state = *st_ptr;
    std::shared_ptr<const State> shst(st_ptr);

    state.t = read_u64();
    state.dimensions = settings.dimensions;
    state.boundary = settings.boundary;

    auto num_readers = read_u32();
    state.readers.reserve(num_readers);
    for (uint32_t i = 0; i < num_readers; i++) {
        state.readers.insert(readReader());
    }

    auto num_books = read_u32();
    state.books.reserve(num_books);
    for (uint32_t i = 0; i < num_books; i++) {
        state.books.insert(readBook());
    }

    return shst;
}

std::pair<eris::eris_id_t, ReaderState> FileStorage::readReader() const {
    auto pair = std::make_pair<eris_id_t, ReaderState>(
            read_u64(),
            ReaderState(settings.dimensions));

    ReaderState &r = pair.second;
    r.id = pair.first;

    // Position
    for (uint32_t d = 0; d < settings.dimensions; d++)
        r.position[d] = read_dbl();

    // Friends
    auto num_friends = read_u32();
    r.friends.reserve(num_friends);
    for (uint32_t i = 0; i < num_friends; i++)
        r.friends.insert(read_u64());

    // Library
    auto libsize = read_u32();
    r.library.reserve(libsize);
    for (uint32_t i = 0; i < libsize; i++) {
        auto status = read_u8();
        auto id = read_u64();
        auto quality = read_dbl();
        r.library.emplace(id, quality);
        bool pirated = status & 1<<1, new_book = status & 1<<2;
        if (status == 1) // Wrote it (other bits not allowed)
            r.wrote.insert(r.wrote.end(), id);
        else if (status & ~(1<<1 | 1<<2))
            throwParseError("found illegal book status (invalid status bits set)");
        else {
            if (new_book) r.new_books.insert(id);
            if (pirated) {
                r.library_pirated.insert(id);
                if (new_book) r.new_pirated.insert(id);
            }
            else {
                r.library_purchased.insert(id);
                if (new_book) r.new_purchased.insert(id);
            }
        }
    }
    // Utility
    r.u = read_dbl();
    r.u_lifetime = read_dbl();
    // Costs
    r.cost_fixed = read_dbl();
    r.cost_unit = read_dbl();
    r.cost_piracy = read_dbl();
    r.income = read_dbl();

    // Beliefs
    belief_data belief = readBelief();
    if (belief.K > 0)
        r.profit = belief.noninformative
            ? Profit(settings.dimensions, belief.K)
            : Profit(settings.dimensions, belief.beta, belief.s2, belief.V, belief.n);

    belief = readBelief();
    if (belief.K > 0)
        r.profit_extrap = belief.noninformative
            ? Profit(settings.dimensions, belief.K)
            : Profit(settings.dimensions, belief.beta, belief.s2, belief.V, belief.n);

    belief = readBelief();
    if (belief.K > 0)
        r.demand = belief.noninformative
            ? Demand(settings.dimensions, belief.K)
            : Demand(settings.dimensions, belief.beta, belief.s2, belief.V, belief.n);

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
    else if (k >= -100 and k < 0) {
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

    // The first K elements are beta values
    belief.beta = VectorXd(k);
    for (unsigned int i = 0; i < belief.K; i++)
        belief.beta[i] = read_dbl();

    // Then s2 and n:
    belief.s2 = read_dbl();
    belief.n = read_dbl();

    // The last K*(K+1)/2 are the V values (but we set them symmetrically in V)
    belief.V = MatrixXd(k, k);
    for (unsigned int r = 0; r < belief.K; r++) {
        for (unsigned int c = r; c < belief.K; c++) {
            double cov = read_dbl();
            belief.V(r,c) = cov;
            belief.V(c,r) = cov;
        }
    }

    return belief;
}

void FileStorage::writeState(const State &state) {
    write_u64(state.t);
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
    write_u64(r.id);

    // Position
    for (size_t i = 0; i < r.position.dimensions; i++)
        write_dbl(r.position[i]);

    // Friends
    write_u32(r.friends.size());
    for (auto &f : r.friends)
        write_u64(f);

    // Library
    write_u32(r.library.size());
    for (auto &l : r.library) {
        uint8_t status = 0;
        if (r.wrote.count(l.first)) status = 1;
        else {
            if (r.library_pirated.count(l.first)) status |= 1<<1;
            if (r.new_books.count(l.first)) status |= 1<<2;
        }
        write_u8(status);
        write_u64(l.first);
        write_dbl(l.second);
    } 
    // Utility
    write_value(r.u);
    write_value(r.u_lifetime);
    // Costs
    write_value(r.cost_fixed);
    write_value(r.cost_unit);
    write_value(r.cost_piracy);
    write_value(r.income);

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

    // Otherwise we're all good, write out the record:
    write_i8(k);

    auto &beta = m.beta();
    // First K elements are the beta values
    for (unsigned int i = 0; i < k; i++)
        write_dbl(beta[i]);

    // Then s2 and n:
    write_dbl(m.s2());
    write_dbl(m.n());

    auto &V = m.V();
    // The last k*(k+1)/2 are the lower triangle of the V matrix, in column major order
    for (unsigned int r = 0; r < k; r++) {
        for (unsigned int c = r; c < k; c++) {
            write_dbl(V(r,c));
        }
    }
}

std::pair<eris_id_t, BookState> FileStorage::readBook() const {
    auto pair = std::make_pair<eris_id_t, BookState>(
            read_u64(), BookState(settings.dimensions));

    BookState &b = pair.second;
    b.id = pair.first;

    b.author = read_u64();
    for (uint32_t d = 0; d < settings.dimensions; d++) {
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
    b.created = read_u64();
    b.lifetime = read_u32();

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
    write_u32(b.sales);
    write_u32(b.sales_lifetime);
    write_u32(b.pirated);
    write_u32(b.pirated_lifetime);
    write_u64(b.created);
    write_u32(b.lifetime);
}

}}
