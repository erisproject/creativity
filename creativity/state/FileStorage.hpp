#pragma once
#include "creativity/state/Storage.hpp"
#include <vector>
#include <fstream>
#include <boost/detail/endian.hpp>

namespace creativity {

class Creativity;

namespace state {

/** Class for file-based storage.  Note that most of the methods of this class will throw errors
 * when the underlying file stream encounters an error.  You should not attempt to use the
 * FileStorage object after such an exception is thrown as the result is unpredictable: the object
 * could simply stop writing new data to the file or the file could become corrupted.
 *
 * Writing new states is done by a background thread which is spawned the first time enqueue() is
 * called with a new state.
 *
 * Note that this object is unable to write simulations with more than 2 billion members (in
 * particular, it stores IDs in the file as 32-bit quantities, and will throw an exception if a
 * value that does not fit in a uint32_t is encountered).
 */
class FileStorage final : public StorageBackend {
    public:
        FileStorage() = delete;

        /** The supported file modes, passed to the constructor. */
        enum class MODE {
            /** Opens the file in read-only mode.  The file must exist and contain valid data,
             * otherwise an exception will be thrown.  Any operation that requires writing to the
             * file (such as adding a State) will fail.
             */
            READONLY,
            /** Opens the file in read-write mode.  If the file exists and is non-empty, it will be
             * parsed (and so must be a valid file).  Otherwise, it will be initialized as a new
             * data file with no states (and uninitialized dimensions and boundaries parameters).
             * New states passed to push_back() will be added to the end of the existing states (if
             * any).
             */
            APPEND,
            /** Just like APPEND, but requires that the file already exist with at least a header.
             */
            READ,
            /** Like APPEND, but truncates any file data first if the file exists (so the file is
             * always treated as an empty file and reinitialized).
             */
            OVERWRITE
        };

        /** Constructs and returns a FileStorage object that uses the given file for reading and
         * (optionally) writing state data.  The file is read or created immediately.
         *
         * \param filename the filename to open
         * \param mode the open mode to use
         *
         * Throws various exceptions if the file does not exist, cannot be read, is empty, or
         * contains invalid data.
         */
        FileStorage(const std::string &filename, MODE mode);

        /** Throws a ParseError exception.  The given message is prefixed with `Parsing file failed
         * [pos=123]: `, where `123` is the current file position.
         */
        void throwParseError(const std::string& message) const;

        /** Exception class thrown when the opened file contains data that cannot be parsed.
         */
        class ParseError : public std::runtime_error {
            public:
                /// Exception constructor
                ParseError(const std::string& message);
        };

        /** Loads the requested state data from the open file into a State object and returns it
         * (wrapped in a shared_ptr).
         */
        virtual std::shared_ptr<const State> load(eris::eris_time_t t) const override;

        /// Returns the number of states currently stored in the file.
        virtual size_t size() const override;

        /** Reads the settings stored in the file header into the given CreativitySettings object.
         */
        virtual void readSettings(CreativitySettings &settings) const override;

        /** Updates the values in the file header to those currently in the CreativitySettings
         * reference provided.
         *
         * Existing values are overwritten.
         */
        virtual void writeSettings(const CreativitySettings &settings) override;

        /** Flushes the file to disc.  This calls `flush()` on the underlying file object.  This is
         * normally not required: any changes will be automatically flushed when the FileStorage
         * object is destroyed; calling this manually primarily guards against unflushed data
         * being unwritten if the process is terminated abnormally (i.e. via a signal or a segfault).
         *
         * Use with caution: frequent file flushing will slow the FileStorage object considerably.
         *
         * \sa std::basic_ostream::flush
         */
        virtual void flush() override;

    protected:

        /// Called from the queue thread to write the given State to the file.
        virtual void thread_insert(std::shared_ptr<const State> &&s) override;

    private:
        /** The file buffer object. Mutable because we need to read from it in const methods. */
        mutable std::fstream f_;

#ifndef CREATIVITY_DISABLE_THREADED_STORAGE
        /** Mutex guarding f_ and related variables (such as state_pos_).  Some operations (such as
         * load/saving settings) are done by the main thread, which accesses f_.
         */
        mutable std::mutex f_mutex_;
#endif

        /// Storage for header data parsing when opening the file, and after writing settings.
        CreativitySettings settings_;

        /** Attempts to parse the metadata (global settings, file locations; essentially everything
         * except for the actual state data) from the opened file.  Throws an exception if parsing
         * fails or the file format appears invalid.
         *
         * This method is called immediately after the file is opened during construction when
         * opened in READONLY or READ modes, and when opened in APPEND mode with a non-empty file.
         */
        void parseMetadata();

        /** Writes a header for an empty FileStorage: that is, a header for a file with 0 states,
         * and with all settings initialized to 0.  After this completes, the file output position
         * will be at the end of the header, i.e. at `FileStorage::HEADER::size`.
         *
         * This method is called during construction when opening a file in OVERWRITE mode or when
         * opening a non-existant or empty file in APPEND mode.
         */
        void writeEmptyHeader();

        /** Reads the given type T from the file at its current location, advancing the file's
         * current location by the size of the given type.
         */
        template <typename T>
        T read_value() const;

        /** Reads a value from the file at its current location into the given variable, advancing
         * the file's current location by the size of the given type.
         */
        template <typename T>
        void read_value(T &val) const { val = read_value<T>(); }

        /** Alias for `read_value<uint64_t>()` */
        uint64_t read_u64() const { return read_value<uint64_t>(); }

        /** Alias for `read_value<uint32_t>()` */
        uint32_t read_u32() const { return read_value<uint32_t>(); }

        /** Alias for `read_value<uint8_t>()` */
        uint8_t read_u8() const { return read_value<uint8_t>(); }

        /** Alias for `read_value<int64_t>()` */
        int64_t read_i64() const { return read_value<int64_t>(); }

        /** Alias for `read_value<int32_t>()` */
        int32_t read_i32() const { return read_value<int32_t>(); }

        /** Alias for `read_value<int8_t>()` */
        int8_t read_i8() const { return read_value<int8_t>(); }

        /** Alias for `read_value<double>()` */
        double read_dbl() const { return read_value<double>(); }

        /** Writes the given value to the file at its current location. */
        template <typename T>
        void write_value(const T &val);

        /** Alias for `write_value((uint64_t) value)` */
        void write_u64(uint64_t value) { write_value(value); }

        /** Alias for `write_value((uint32_t) value)` */
        void write_u32(uint32_t value) { write_value(value); }

        /** Alias for `write_value((uint8_t) value)` */
        void write_u8(uint8_t value) { write_value(value); }

        /** Alias for `write_value((int64_t) value)` */
        void write_i64(int64_t value) { write_value(value); }

        /** Alias for `write_value((int32_t) value)` */
        void write_i32(int32_t value) { write_value(value); }

        /** Alias for `write_value((int8_t) value)` */
        void write_i8(int8_t value) { write_value(value); }

        /** Alias for `write_value((double) value)` */
        void write_dbl(double value) { write_value(value); }

        /** Seeks to the end of the file, where a new block can be written.  If required, padding is
         * added to the end of the file so that the new block is on a block alignment boundary.  The
         * file output position will be at the end of the file after this call; the location is also
         * returned for convenience.
         */
        int64_t newBlock();

        /** Copies the appropriate number of bytes from memory from the location of `from` into a
         * variable of type T and returns it.  This is safer than a reinterpret_cast as it avoids
         * alignment issues (i.e. where `from` has a memory alignment that is invalid for T
         * variables).
         *
         * The number of bytes copied is `sizeof(T)`; `from` may be a smaller type (such as a
         * `char`) as long as `sizeof(T)` bytes of valid data exists at its location.
         */
        template <typename T, typename F>
        static T parse_value(const F &from);

        /// Like parse_value(from), but copies into the given variable instead of returning it.
        template <typename T, typename F>
        static void parse_value(const F &from, T &to) { to = parse_value<T>(from); }

        /** Copies the appropriate number of bytes from memory from the location of `from` into the
         * location of `to`, reordered (if necessary) to be in the file's little-endian order.
         *
         * The number of bytes copied is `sizeof(F)`; `to` may be a smaller type (such as a `char`)
         * as long as `sizeof(F)` bytes of valid space exists at its location.
         */
        template <typename F, typename T>
        static void store_value(const F &from, T &to);

        /// Store the locations of the state data for each state;
        std::vector<std::streampos> state_pos_;

        /** The header contains the first 485 state locations (bytes 208 through 4887); if the file
         * has more states, the location at 4088 points to a "continuation block": 511 state
         * locations and a further continuation block location.  This variable stores the file
         * locations of continuation blocks, read when the file is first opened and updated as
         * updates require additional blocks.
         */
        std::vector<std::streampos> cont_pos_;

        /** Everything (header, library blocks, continuation blocks, states) gets aligned to this size,
         * with padding added as required before beginning a new block.
         */
        static constexpr unsigned int BLOCK_SIZE = 4096;

        /// Constants for header attributes
        struct HEADER {
            static constexpr unsigned int size = BLOCK_SIZE; ///< Size of the header, in bytes
            /// The magic value 'CrSt' that identifies the file as created by this class
            static constexpr char fileid[4] = {'C', 'r', 'S', 't'};
            /// The test value bytes, which should be interpreted as the various test values below
            static constexpr char test_value[8] = {20,106,10,-50,0,24,69,-64};
            /// Header positions of various fields
            struct pos {
                static constexpr int64_t
                    fileid = 0, ///< First 4 bytes are File ID such as "CrSt" ("Cr"eative "St"ate file)
                    filever = 4, ///< Second 4 bytes are u32 file format version
                    test_value = 8, ///< the 8-byte test value (which is interpreted as various types)
                    num_states = 16, ///< the number of states stored in this file (u32)
                    dimensions = 20, ///< the number of dimensions of the simulation these states belong to (u32)
                    readers = 24, ///< the number of readers in the simulation
                    boundary = 28, ///< the simulation boundary (positive double)
                    book_distance_sd = 36, ///< standard deviation of distance of new books from authors (positive double)
                    book_quality_sd = 44, ///< standard deviation of perceived book quality draw (draw is normal, with mean = true quality) (dbl)
                    reader_step_sd = 52, ///< reader random step sd (dbl)
                    reader_creation_shape = 60, ///< reader q(l) shape parameter (dbl)
                    reader_creation_scale_min = 68, ///< reader q(l) scale parameter ~ U[a,b]; this is 'a' (dbl)
                    reader_creation_scale_max = 76, ///< reader q(l) scale parameter ~ U[a,b]; this is 'b' (dbl)
                    creation_time = 84, ///< Time to create a book
                    cost_fixed = 88, ///< Fixed cost of keeping a book on the market (dbl)
                    cost_unit = 96, ///< Unit cost of an author creating a copy of a book (dbl)
                    cost_piracy = 104, ///< Unit cost of getting a copy of a book via piracy (dbl)
                    income = 112, ///< Per-period reader external income (before incurring authorship costs or receiving book profits) (dbl)
                    piracy_begins = 120, ///< the sharing start period (eris_time_t = u32)
                    piracy_link_proportion = 124, ///< The proportion of potential friendship links that exist
                    prior_scale = 132, ///< The prior 'n' multiplier (usually 0-1)
                    prior_scale_piracy = 140, ///< The prior 'n' multiplier in the first piracy period
                    prior_scale_burnin = 148, ///< The prior 'n' multiplier (usually 0-1)
                    burnin_periods = 156, ///< The burnin-periods during which beliefs are more heavily discounted (u32)
                    init_prob_write = 160, ///< The probability of writing (while beliefs noninformative)
                    init_q_min = 168, ///< `a` in U[a,b], the noninformative belief authorship quality level
                    init_q_max = 176, ///< `b` in U[a,b], the noninformative belief authorship quality level
                    init_p_min = 184, ///< `a` in U[a,b], the noninformative belief book price
                    init_p_max = 192, ///< `b` in U[a,b], the noninformative belief book price
                    init_prob_keep = 200, ///< The probability of keeping a book on the market (uninformed beliefs)
                    init_keep_price = 208, ///< If keeping a book on the market, the new price is (p-c)*s+c, where this is s
                    init_belief_threshold = 216, ///< The required n-k value for readers to use beliefs instead of initial behaviour (i32)
                    // NB: this next one should be an integer multiple of 8 (so that states+cont
                    // location go to the end, and so the `states` division below has no remainder)
                    //
                    // padding, if needed, here.
                    state_first = 224, ///< the first state record
                    state_last = size - 16, ///< the last state record
                    continuation = size - 8; ///< the header continuation block pointer (used once header state blocks fill up)
            };
            /** The number of states that can be stored in the header.  Storing additional states
             * requires using continuation blocks.
             */
            static constexpr unsigned int states = (pos::state_last - pos::state_first) / 8 + 1;
            static constexpr uint32_t u32_test = 3456789012; ///< unsigned 32-bit integer test value
            static constexpr int32_t  i32_test = -838178284; ///< signed 32-bit integer test value
            static constexpr uint64_t u64_test = 13854506220411054612ul; ///< unsigned 64-bit integer test value
            static constexpr int64_t  i64_test = -4592237853298497004; ///< signed 64-bit integer test value
            // The following constant is exactly representable in a 64-bit double (in fact it is the
            // exact value of the test_value characters above, interpreted as a little-endian stored
            // IEEE754 double).
            static constexpr double dbl_test = -42.187524561963215319337905384600162506103515625; ///< double test value
        };

        /// Constants for continuation blocks
        struct CBLOCK {
            static constexpr unsigned int states = (BLOCK_SIZE / 8) - 1; ///< Number of states stored in a continuation block
            static constexpr unsigned int next_cblock = states * 8; ///< Relative position of the pointer to the next continuation block
        };

        /// A BLOCK_SIZE block of 0s
        static constexpr char ZERO_BLOCK[BLOCK_SIZE] = {0};

        /** Constants for library storage.  Libraries are stored in blocks of many packed
         * (bookid,t,quality,status) tuples followed by a pointer to the next library block.  When
         * reading a reader state, the first block of the library is pointed to, and is read until
         * either a 0 id or a `t` value larger than that of the state being read is found.
         */
        struct LIBRARY {
            /** a reader library record size: the book id (u32), the acquisition period (u32), the
             * reader-specific quality (dbl), and the library book status (u8: 0=wrote, 1=bought,
             * 2=pirated).
             */
            static constexpr unsigned int record_size = sizeof(uint32_t) + sizeof(uint32_t) + sizeof(double) + sizeof(uint8_t);
            static constexpr unsigned int block_records = (BLOCK_SIZE - 8)/record_size; ///< Number of records that will fit in a block
        };

        /// Flags for opening in read-only mode
        static constexpr std::ios_base::openmode open_readonly = std::ios_base::binary | std::ios_base::in;
        /// Flags for opening in overwrite mode
        static constexpr std::ios_base::openmode open_overwrite = open_readonly | std::ios_base::out | std::ios_base::trunc;
        /** Flags for opening for read-write when the file exists (if it doesn't exist,
         * open_overwrite is used instead because this mode fails when the file doesn't exist).
         */
        static constexpr std::ios_base::openmode open_readwrite = open_readonly | std::ios_base::out;

        /** Parses the requested number of std::streampos values from memory beginning at `from` and
         * adds them to `state_pos_`.  Throws an exception if the read fails or values are invalid
         * (less than `FileStorage::HEADER::size` >= `end`).  Note that this simply reads a
         * contiguous set of state locations, it does *not* handle continuation block logic (but
         * rather continuation block reading uses this method).
         *
         * \param from reference to the variable located at the first state position value
         * \param count the number of state locations to read
         * \param end the maximum+1 state location that will be accepted; typically the current file
         * size
         */
        void parseStateLocations(const char &from, const size_t count, const size_t end);

        /** Adds the given location to the header or last continuation block, as appropriate.
         * Creates a new continuation block (updating the previous one) if necessary.
         *
         * After updating (or creating) the continuation block(s), this also updates the header with
         * the new number of states, and adds the location to state_pos_.
         */
        void addStateLocation(std::streampos location);

        /** Adds a continuation block to the end of the file and updates the last current
         * continuation block (or, if none, the header) to point to this new continuation block.
         * This also adds an element to cont_pos_.
         *
         * This method should only be called when the current last continuation block (or the states
         * in the header) is full.
         *
         * The file output position is not guaranteed to be at any particular position after this
         * call.
         */
        void createContinuationBlock();

        // Sorting method that sorts by the library item's acquired date
        class lib_comp_less final {
            public:
                using ref = const std::pair<const uint32_t, BookCopy>&;
                bool operator()(ref a, ref b) const {
                    return a.second.acquired < b.second.acquired;
                }
        };

        // Reader library data
        struct lib_data {
            int64_t pos_next; // Next library row position
            unsigned int records_remaining; // How many library rows are remaining (0 = need a new block)
            // book id to BookCopy
            std::unordered_map<uint32_t, BookCopy> library;
            // library, but sorted by acquired date (for fast retrieval)
            std::multiset<std::reference_wrapper<std::pair<const uint32_t, BookCopy>>, lib_comp_less> library_acq_sorted;
            lib_data(int64_t pos_next) : pos_next(pos_next), records_remaining(LIBRARY::block_records) {}
        };

        /** Reader library data.  Key is the reader id, value is the data. */
        std::unordered_map<unsigned int, lib_data> reader_lib_;

        /** Read the reader library pointer block, which immediately follows the header.
         *
         * The block consists of num_reader (reader id, library block pointer) tuples.
         *
         * The pointed at library block consists of (bookid,t,quality,status) tuples for all books
         * in the reader's library over all stored simulation periods.  If the block overflows,
         * additional continued blocks are created and pointed at by a pointer at the end of the
         * library block.
         *
         * Reader blocks are stored in memory when the file is opened and updated (both in memory
         * and on disk) as states are added.
         */
        void parseLibraryPointerBlock();

        /** Writes the reader library pointer block.  This is called when the very first state is
         * added to a new file, and as a result will end up putting the reader pointer block
         * immediately after the header.  The keys of the given map are used to write the locations.
         */
        void writeLibraryPointerBlock(const std::unordered_map<eris::eris_id_t, ReaderState> &readers);

        /** Ensures that the library block for reader `r` has all of the reader's library books in
         * it.  If any are missing, they are added to disk and to reader_lib_.
         */
        void updateLibrary(const ReaderState &r);

        /** Creates an empty library block at the end of the file and returns the location of the
         * beginning of the block.
         *
         * The file output position is not guaranteed to be at any particular position after this
         * call.
         */
        int64_t createLibraryBlock();

        /** Reads a state starting at the current file position.  The state data structure is as
         * follows:
         *
         *     u32      t (simulation time period)
         *     READER[] readers array
         *     BOOK[]   books
         *
         * where type[] indicates an array structured as:
         *     u32          length
         *     type*length  sequential type records
         *
         * \see readReader() for the structure of READER
         * \see readBook() for the structure of BOOK
         */
        std::shared_ptr<const State> readState() const;

        /// Writes the given state at the current file position.
        void writeState(const State &state);

        /** Reads a ReaderState record from the current file position and returns it in an
         * {eris_id_t, ReaderState} pair, where `.first` is the id.  Such a record consists of:
         *
         *     u32              id
         *     dbl*DIM          position (DIM = dimensions)
         *     u32[]            friend ids
         *     dbl              u
         *     dbl              u_lifetime
         *     dbl              cost_fixed
         *     dbl              cost_unit
         *     dbl              cost_piracy
         *     dbl              income
         *     dbl              creation_shape
         *     dbl              creation_scale
         *     BELIEF           profit belief
         *     BELIEF           profit extrapolated belief
         *     BELIEF           demand belief
         *     BELIEF           quality belief
         *     BELIEF[]         profit stream beliefs (for different K() values)
         *
         * where type[] indicates an array structured as:
         *     u32          length
         *     type*length  sequential type records
         *
         * and (type1,type2) indicates a single type1 value followed immediately by a single type2
         * value.
         *
         * The reader's library is also read, but is not contained in the reader block; rather it is
         * stored in the locations referenced in the dedicated reader library section at the
         * beginning of the file.
         *
         * BELIEF is a set of belief data, as handled by readBelief(int64_t).  profit extrapolated
         * belief is set to a belief only if it differs from profit_belief; otherwise it is simply
         * set to a no-data, noninformative belief record (and shouldn't be used).
         *
         * profit stream beliefs may not be placeholder beliefs (i.e. default constructed objects);
         * such objects should simply be omitted when writing the data.  Each profit stream belief K
         * value must also be unique.
         */
        std::pair<eris::eris_id_t, ReaderState> readReader(eris::eris_time_t t) const;

        /// Writes a reader state at the current file position.
        void writeReader(const ReaderState& reader);

        /// Structure holding parsed belief data
        typedef struct {
            uint32_t K; ///< Number of parameters; K=0 for a default constructed (invalid) model (the remaining values will be uninitialized)
            bool noninformative; ///< True if this is a noninformative model (in which case the following are not set)
            Eigen::VectorXd beta; ///< eris::belief::BayesianLinear beta vector
            double s2; ///< eris::belief::BayesianLinear s2 value
            double n; ///< eris::belief::BayesianLinear n value
            Eigen::MatrixXd Vinv; ///< eris::belief::BayesianLinear Vinv matrix
            bool draw_gibbs; ///< If true, the last draw from this belief used Gibbs sampling (false = no draws, or rejection sampling)
            uint32_t draw_success_cumulative, ///< For a restricted belief, the number of successful draws
                     draw_discards_cumulative; ///< For a restricted belief, the number of discarded draws
        } belief_data;

        /** Reads a belief structure from the current file location.  The structure is as follows:
         *
         * In the case of a default-constructed model or a noninformative model, no further reading
         * is required.  When a file location (or the "immediately following" -512 value) is
         * recognized, that indicated location contains a belief record structured as follows:
         *
         *     i8       K, the number of model parameters
         *     u8       status bits (described below)
         * then (if noninformative bit not set):
         *     dbl*K    beta vector (K values)
         *     dbl      s2
         *     dbl      n
         *     dbl*Z    the `Z=K(K+1)/2` elements of the lower triangle of the V^{-1} matrix, in column-major
         *              order (that is, [0,0], [1,0], [2,0], ..., [K,0], [1,1], [2,1], ..., [K,1], ..., [K,K])
         *     u32      the cumulative successful LinearRestricted draws (only for LinearRestricted
         *              models)
         *     u32      the cumulative discarded LinearRestricted draws (only for LinearRestricted
         *              models)
         *
         * The K value has the following interpretations: if in [1,120], the model has this number of
         * parameters, and is not a completely noninformative model.  If -128, the model object is
         * default-constructed (and thus has no useful data).  Anything else is invalid.
         *
         * The second byte is a bit field, which currently has the following interpretations:
         * - lowest bit (bit & 1): the model is a LinearRestricted model, and thus carries extra
         *   LinearRestricted data (the two u32's above).  If not set, the u32s are not present in
         *   the record.
         * - bit 2 (bit & 2): the last draw from this model used Gibbs sampling (only applicable to
         *   restricted models)
         * - bit 3 (bit & 4): the model is noninformative.  If this bit is set, then instead of
         *   beta/V/etc. data the record contains non-informative data that hasn't been loaded yet;
         *   otherwise the beta/V/n/s2 data follows.
         * - other bits are currently unused.
         */
        belief_data readBelief() const;

        /** Writes a belief at the current file location.  The generally consists of a single i64
         * value at the current file location, and, if required, a belief record immediately
         * following it.
         *
         * If the belief exactly matches a belief that has been written or read by the current
         * object, this reuses that record (by just writing the existing location).  This may still
         * result in duplicate records if the file was loaded from an existing state file as state
         * data is only loaded as needed: deduplication can only occur for known beliefs.
         *
         * \sa readBelief() for the structure and values actually written
         */
        void writeBelief(const eris::belief::BayesianLinear &belief);

        /** Reads a book state data from the current file position and returns it in a <eris_id_t,
         * BookState> pair, where the eris_id_t is the book id.  The book data is:
         *
         *     u32          id
         *     u32          author id
         *     dbl*DIM      position (DIM = dimensions)
         *     dbl          quality
         *     dbl          price (market is derived from this: market=true unless price is NaN)
         *     dbl          revenue
         *     dbl          revenueLifetime
         *     u32          sales
         *     u32          salesLifetime
         *     u32          pirated
         *     u32          piratedLifetime
         *     u32          created
         *     u32          lifetime
         */
        std::pair<eris::eris_id_t, BookState> readBook() const;

        /** Writes a book at the current file position.  See readBook() for data layout. */
        void writeBook(const BookState &book);

};

#if defined(BOOST_BIG_ENDIAN)
#include "creativity/state/FileStorage-BE.hpp"
#elif defined(BOOST_LITTLE_ENDIAN)
#include "creativity/state/FileStorage-LE.hpp"
#else
#error System endianness not supported (neither big-byte nor little-byte endianness detected)!
#endif

}}
