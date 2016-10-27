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
using namespace eris::learning;
using namespace creativity::belief;
using namespace Eigen;
using DrawMode = BayesianLinearRestricted::DrawMode;
namespace fs = boost::filesystem;

FileStorage::FileStorage(CreativitySettings &settings, std::unique_ptr<std::stringstream> &&s, Mode mode)
    : StorageBackend(settings)
{
    memory(std::move(s), mode);
}


void FileStorage::throwParseError(const std::string& message) const {
    decltype(f_->tellg()) pos;
    try { pos = f_->tellg(); }
    catch (std::ios_base::failure &f) { pos = -1; } // -1 is also be returned by .tellg() for pre-existing failures
    throw parse_error("Parsing file failed [pos=" + (pos == -1 ? std::string{"(error)"} : std::to_string(pos)) + "]: " + message);
}

void FileStorage::readExtraHeader() {
    state_pointer_block_ = f_->tellg();
    blockListSkip(state_pointer_block_);
    readerlib_block_ = f_->tellg();
    blockListIterate(readerlib_block_, [&] (uint32_t) -> bool {
            auto readerid = read<eris_id_t>();
            auto &rlib = reader_lib_[readerid];
            rlib.first = read<uint64_t>();
            // Iterate through the readers library and store everything in it:
            blockListIterate(rlib.first, [&] (uint32_t) -> bool {
                    // NB: if this changes, remember to update updateLibraries as well.
                    auto book_id = read<eris_id_t>();
                    auto quality = read<double>();
                    auto status_code = read<uint8_t>();
                    auto acquired = read<eris_time_t>();

                    if (status_code > 3) throw std::runtime_error("Found unknown BookCopy status code `" + std::to_string(status_code) + "'");
                    BookCopy::Status status = 
                        status_code == 0 ? creativity::BookCopy::Status::wrote :
                        status_code == 1 ? creativity::BookCopy::Status::purchased_market :
                        status_code == 2 ? creativity::BookCopy::Status::pirated :
                        creativity::BookCopy::Status::purchased_public;
                    auto ins = rlib.second.library.emplace(std::move(book_id), BookCopy(quality, status, acquired));
                    rlib.second.library_acq_sorted.emplace(std::ref(*ins.first));
                    return true;
            });

            return true;
    });
}

void FileStorage::writeExtraHeader() {
    state_pointer_block_ = pointerListCreate();
    constexpr auto size = sizeof(eris_id_t) + sizeof(uint64_t);
    readerlib_block_ = blockListCreate(size, 4088/size); // Fit into 4KiB (including the pointer)
}

size_t FileStorage::size() {
    return pointerListSize(state_pointer_block_);
}

void FileStorage::thread_insert(std::shared_ptr<const State> &&state) {
    // First we need to update the libraries of readers, adding any newly acquired items:
    updateLibraries(state->readers);

    // Ensure that the state at index t is really state t
    auto sz = size();
    if (sz != state->t) throw std::runtime_error("States missing/out of order: trying to add state " + std::to_string(state->t) + ", but file only contains " + std::to_string(sz) + " states");

    // Now seek to the end of the file to write our state record.
    f_->seekp(0, f_->end);
    uint64_t location = f_->tellp();
    writeState(*state);
    pointerListAppend(state_pointer_block_, location);
}


std::shared_ptr<const State> FileStorage::load(eris_time_t t) {
    try {
        pointerListSeek(state_pointer_block_, t);
    }
    catch (std::out_of_range&) {
        return std::shared_ptr<const State>();
    }

    return readState();
}

void FileStorage::storage_flush() {
    // Ensure the filehandle output is flushed to disk
    if (f_) f_->flush();
}

void FileStorage::configureHeaderFields() {
    addHeaderField(settings_.dimensions);
    addHeaderField(settings_.readers);
    addHeaderField(settings_.boundary);
    addHeaderField(settings_.book_distance_mean);
    addHeaderField(settings_.book_quality_sd);
    addHeaderField(settings_.reader_step_mean);
    addHeaderField(settings_.reader_creation_shape);
    addHeaderField(settings_.reader_creation_scale_min);
    addHeaderField(settings_.reader_creation_scale_range);
    addHeaderField(settings_.creation_time);
    addHeaderField(settings_.creation_fixed);
    addHeaderField(settings_.cost_market);
    addHeaderField(settings_.cost_unit);
    addHeaderField(settings_.cost_piracy);
    addHeaderField(settings_.income);
    addHeaderField(settings_.piracy_begins);
    addHeaderField(settings_.piracy_link_proportion);

    addHeaderField(settings_.policy);
    addHeaderField(settings_.policy_begins);

    addHeaderField(settings_.policy_public_sharing_tax);

    addHeaderField(settings_.policy_public_sharing_voting_tax);
    addHeaderField(settings_.policy_public_sharing_voting_votes);

    addHeaderField(settings_.policy_catch_tax);
    addHeaderField(settings_.policy_catch_cost);
    addHeaderField(settings_.policy_catch_fine);
    addHeaderField(settings_.policy_catch_mu);
    addHeaderField(settings_.policy_catch_sigma);

    addHeaderField(settings_.prior_scale);
    addHeaderField(settings_.prior_scale_burnin);
    addHeaderField(settings_.prior_scale_piracy);
    addHeaderField(settings_.burnin_periods);
    addHeaderField(settings_.initial.prob_write);
    addHeaderField(settings_.initial.l_min);
    addHeaderField(settings_.initial.l_range);
    addHeaderField(settings_.initial.p_min);
    addHeaderField(settings_.initial.p_range);
    addHeaderField(settings_.initial.prob_keep);
    addHeaderField(settings_.initial.keep_price);
    addHeaderField(settings_.initial.belief_threshold);
}

void FileStorage::writeSettings() {
    updateHeaderFields();
}

void FileStorage::updateLibraries(const std::map<eris_id_t, ReaderState> &readers) {
    for (const auto &rp : readers) {
        const eris_id_t &rid = rp.first;
        const auto &r = rp.second;
        auto lib = reader_lib_.find(rid);
        uint64_t lib_loc = 0;

        if (lib == reader_lib_.end()) {
            // The reader wasn't found, so create a library block for it
            constexpr auto libsize = sizeof(eris_id_t) + sizeof(double) + sizeof(uint8_t) + sizeof(eris_time_t);
            uint64_t lib_loc = blockListCreate(libsize, 4088/libsize); // Aim for ~4KiB library blocks (including pointer)

            // Now add the reader and library pointer:
            blockListAppend(readerlib_block_);
            *this << rid << lib_loc;

            auto &p = reader_lib_[rid];
            p.first = lib_loc;
        }
        else {
            lib_loc = lib->second.first;
        }

        // Now copy any new/missing library values into the library block
        lib_data &lib_stored = reader_lib_[rid].second;
        for (auto &l : r.library) {
            if (not lib_stored.library.count(l.first)) {
                blockListAppend(lib_loc);
                auto &bc = l.second;
                *this << l.first
                    << bc.quality
                    << uint8_t(bc.wrote() ? 0 : bc.purchased_market() ? 1 : bc.pirated() ? 2 : /* bc.purchased_public() */ 3)
                    << bc.acquired;

                auto ins = lib_stored.library.emplace(l.first, l.second);
                lib_stored.library_acq_sorted.emplace(std::ref(*ins.first));
            }
        }
    }
}

std::shared_ptr<const State> FileStorage::readState() {
    State *st_ptr = new State();
    State &state = *st_ptr;
    std::shared_ptr<const State> shst(st_ptr);

    *this >> state.t;
    state.dimensions = settings_.dimensions;
    state.boundary = settings_.boundary;

    for (auto type = read<uint8_t>(); type != TYPE_DONE; *this >> type) {
        switch (type) {
            case TYPE_READERS:
                {
                    auto num = read<uint32_t>();
                    for (uint32_t i = 0; i < num; i++) state.readers.insert(readReader(state.t));
                    break;
                }
            case TYPE_BOOKS:
                {
                    auto num = read<uint32_t>();
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

void FileStorage::writeState(const State &state) {
    if (state.readers.size() != settings_.readers)
        throw std::runtime_error("FileStorage error: cannot write a simulation where the number of readers changes");

    // First the state's `t`:
    *this << state.t;

    // NB: The order of these doesn't technically matter, but state.readers is an ordered map so
    // that we order by eris_id_t and thus get reproducible files for identical state data.
    if (not state.readers.empty()) {
        *this << uint8_t(TYPE_READERS) << uint32_t(state.readers.size());
        for (auto &r : state.readers) writeReader(r.second);
    }

    if (not state.books.empty()) {
        *this << uint8_t(TYPE_BOOKS) << uint32_t(state.books.size());
        for (auto &b : state.books) writeBook(b.second);
    }

    if (state.publictracker and state.publictracker->id != 0) {
        *this << uint8_t(TYPE_PUBLIC_TRACKER);
        writePublicTracker(*state.publictracker);
    }

   *this << uint8_t(TYPE_DONE);
}


std::pair<eris_id_t, ReaderState> FileStorage::readReader(eris_time_t t) {
    auto pair = std::make_pair<eris_id_t, ReaderState>(
            read<eris_id_t>(),
            ReaderState(settings_.dimensions));

    ReaderState &r = pair.second;
    r.id = pair.first;

    // Position
    for (unsigned d = 0; d < settings_.dimensions; d++)
        *this >> r.position[d];

    // Friends
    auto num_friends = read<uint32_t>();
    for (uint32_t i = 0; i < num_friends; i++)
        r.friends.emplace(read<eris_id_t>());

    // Library (sorted by acquisition period):
    auto &lib_sorted = reader_lib_.at(r.id).second.library_acq_sorted;
    for (auto &l : lib_sorted) {
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

    // Utility and creation parameters:
    *this >> r.u >> r.u_lifetime;
    *this >> r.creation_shape >> r.creation_scale;

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

    auto pstream_locs = read<uint32_t>();
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

void FileStorage::writeReader(const ReaderState &r) {
    FILESTORAGE_DEBUG_WRITE_START
    *this << r.id;

    // Position
    for (size_t i = 0; i < r.position.dimensions; i++)
        *this << r.position[i];

    // Friends
    *this << uint32_t(r.friends.size());
    for (auto &f : r.friends)
        *this << f;

    // Utility
    *this << r.u << r.u_lifetime;
    // Creation parameters
    *this << r.creation_shape << r.creation_scale;

    FILESTORAGE_DEBUG_WRITE_CHECK(8+8*settings_.dimensions+4+8*r.friends.size()+8+8+8+8)

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
    *this << ps_count;
    for (auto &b : r.profit_stream) {
        if (b.second.K() > 0) writeBelief(b.second);
    }
}


FileStorage::belief_data FileStorage::readBelief() {
    belief_data belief;
    auto k = read<int8_t>();
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
    auto status = read<uint8_t>();
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
            *this >> belief.beta[i];

        // Then s2 and n:
        *this >> belief.s2 >> belief.n;

        // Then K*(K+1)/2 V values (but we set them symmetrically in V)
        belief.Vinv = MatrixXd(k, k);
        for (unsigned int r = 0; r < belief.K; r++) {
            for (unsigned int c = 0; c <= r; c++) {
                double cov = read<double>();
                belief.Vinv(r,c) = cov;
                if (c != r) belief.Vinv(c,r) = cov;
            }
        }

        // If this was a restricted model, we read the draw discards and success values
        if (restricted_model) {
            belief.draw_success_cumulative = read<uint32_t>();
            belief.draw_discards_cumulative = read<uint32_t>();
        }
    }

    return belief;
}

void FileStorage::writeBelief(const BayesianLinear &m) {
    FILESTORAGE_DEBUG_WRITE_START
    auto &k = m.K();
    if (k > 120) {
        throw std::runtime_error("creativity::state::FileStorage cannot handle beliefs with K > 120");
    }
    else if (k == 0) {
        // If K = 0 (a default-constructed, non-model object), we just write out the special value and we're done
        *this << int8_t(-128);
        return;
    }

    // Otherwise write out the full record, starting with the positive size:
    *this << int8_t(k);

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
    *this << status;

    if (not m.noninformative()) {
        auto &beta = m.beta();
        // First K elements are the beta values
        for (unsigned int i = 0; i < k; i++)
            *this << beta[i];

        // Then s2 and n:
        *this << m.s2() << m.n();

        const auto &Vinv = m.Vinv();
        // The last k*(k+1)/2 are the lower triangle of the V matrix, in row major order
        for (unsigned int r = 0; r < k; r++) {
            for (unsigned int c = 0; c <= r; c++) {
                *this << Vinv(r,c);
            }
        }

        if (restricted_model) {
            *this << uint32_t(lr->draw_rejection_success) << uint32_t(lr->draw_rejection_discards);
        }
    }

    FILESTORAGE_DEBUG_WRITE_CHECK(
            2 + (
                m.noninformative()
                ? 0 : (8*(2+k+k*(k+1)/2) + (restricted_model ? 4*2 : 0))
                )
            );
}

std::pair<eris_id_t, BookState> FileStorage::readBook() {
    auto pair = std::make_pair<eris_id_t, BookState>(
            read<eris_id_t>(), BookState(settings_.dimensions));

    BookState &b = pair.second;
    b.id = pair.first;

    *this >> b.author;
    for (uint32_t d = 0; d < settings_.dimensions; d++) {
        *this >> b.position[d];
    }

    *this >> b.quality;
    auto stat = read<uint8_t>();
    b.market_private = stat & 1;
    stat &= ~1;
    if (stat != 0) throw std::runtime_error("FileStorage error: found invalid/unknown book status bit");
    *this >> b.price >> b.revenue >> b.revenue_lifetime >> b.prize >> b.prize_lifetime
        >> b.sales >> b.sales_lifetime_private >> b.sales_lifetime_public >> b.pirated >> b.pirated_lifetime
        >> b.created >> b.lifetime_private >> b.votes >> b.votes_lifetime;

    return pair;
}

void FileStorage::writeBook(const BookState &b) {
    FILESTORAGE_DEBUG_WRITE_START

    *this << b.id << b.author;

    for (size_t i = 0; i < b.position.dimensions; i++)
        *this << b.position[i];

    *this << b.quality;
    uint8_t status = 0;
    if (b.market_private) status |= 1;
    *this << status;
    *this << b.price << b.revenue << b.revenue_lifetime << b.prize << b.prize_lifetime
        << b.sales << b.sales_lifetime_private << b.sales_lifetime_public << b.pirated << b.pirated_lifetime
        << b.created << b.lifetime_private << b.votes << b.votes_lifetime;

    FILESTORAGE_DEBUG_WRITE_CHECK(8+8+8*settings_.dimensions+8+1+8*5+4*9);
}

std::unique_ptr<PublicTrackerState> FileStorage::readPublicTracker() {
    std::unique_ptr<PublicTrackerState> ptsptr(new PublicTrackerState);
    auto &pts = *ptsptr;
    *this >> pts.id >> pts.dl_tax >> pts.vote_tax >> pts.dl_unspent >> pts.vote_unspent;
    return ptsptr;
}

void FileStorage::writePublicTracker(const PublicTrackerState &pts) {
    *this << pts.id << pts.dl_tax << pts.vote_tax << pts.dl_unspent << pts.vote_unspent;
}


}}
