#pragma once
#include "creativity/state/Storage.hpp"
#include <vector>
#include <fstream>
#include <boost/predef/other/endian.h>

namespace creativity { namespace state {

/** Class for file-based storage.  Note that most of the methods of this class will throw errors
 * when the underlying file stream encounters an error.  You should not attempt to use the
 * FileStorage object after such an exception is thrown as the result is unpredictable: the object
 * could simply stop writing new data to the file or the file could become corrupted.
 *
 * This object is also not thread-safe on its own: you must ensure that only one thread at a time
 * accesses it, even for read-only methods.
 */
class FileStorage : public Storage, private eris::noncopyable {
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
         * (optionally) writing state data.
         *
         * \param filename the filename to open
         * \param mode the open mode to use
         *
         * Throws various exceptions if the file does not exist, cannot be read, is empty, or
         * contains invalid data.
         */
        FileStorage(const std::string &filename, MODE mode);

        /** Throws a ParseError exception.  The given message is prefixed with `Parsing file failed
         * [pos=123]: `, where `123` is the current f_ position.
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
        virtual std::shared_ptr<const State> operator[](size_t i) const override;

        /// Returns the number of states currently stored in the file.
        virtual size_t size() const override;

        /// Adds a new State to the file.
        virtual void push_back(std::shared_ptr<const State> state) override;

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
        /** The file buffer object. Mutable because we need to read from it in const methods. */
        mutable std::fstream f_;

        /** Attempts to parse the metadata (global settings, file locations; essentially everything
         * except for the actual state data) from the opened file.  Throws an exception if parsing
         * fails or the file format appears invalid.
         *
         * This method is called immediately after the file is opened during construction when
         * opened in READONLY or READ modes, and when opened in APPEND mode with a non-empty file.
         */
        void parseMetadata();

        /** Writes a header for an empty FileStorage: that is, a header for a file with 0 states.
         * After this completes, the file output position will be at the end of the header, i.e. at
         * `FileStorage::HEADER::size`.
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

        /** Alias for `write_value((double) value)` */
        void write_dbl(double value) { write_value(value); }

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

        /** Copies the appropriate number of bytes from memory from the location of `from` into the
         * location of `to`, reordered (if necessary) to be in the file's little-endian order.
         *
         * The number of bytes copied is `sizeof(F)`; `to` may be a smaller type (such as a `char`)
         * as long as `sizeof(T)` bytes of valid space exists at its location.
         */
        template <typename F, typename T>
        static void store_value(const F &from, T &to);

        /// Store the locations of the state data for each state;
        std::vector<std::streampos> state_pos_;

        /** Stores weak references of parsed State data; as long as a returned State is still
         * referenced somewhere else, we can simply return a new reference to it; otherwise it needs
         * to be recreated.
         */
        mutable std::vector<std::weak_ptr<const State>> states_;

        /** The header contains the first 57 state locations (bytes 48 through 504); if the file has
         * more states, the location at 504 points to a "continuation block": 63 state locations and
         * a further continuation block location.  This variable stores the file locations of
         * continuation blocks, read when the file is first opened and updated as updates require
         * additional blocks.
         */
        std::vector<std::streampos> cont_pos_;

        /// Constants for header attributes
        struct HEADER {
            static constexpr unsigned int size = 512; ///< Size of the header, in bytes
            /// The magic value 'CrSt' that identifies the file as created by this class
            static constexpr char fileid[4] = {'C', 'r', 'S', 't'};
            /// The test value bytes, which should be interpreted as the various test values below
            static constexpr char test_value[8] = {20,106,10,-50,0,24,69,-64};
            /// Header positions of various fields
            struct pos {
                static constexpr int64_t
                    fileid = 0, ///< First 4 bytes are File ID such as "CrSt" ("Cr"eative "St"ate file)
                    filever = 4, ///< Second 4 bytes are u32 file format version
                    test_value = 8, ///< location of the 8-byte test value (which is interpreted as various types)
                    states = 16, ///< location of the number of states stored in this file (u32)
                    dimensions = 20, ///< number of dimensions of the simulation these states belong to (u32)
                    boundary = 24, ///< Location of the simulation boundary (positive double)
                    sharing_begins = 32, ///< Location of the sharing start period (u64)
                    state_first = 40, ///< Location of the first state record
                    state_last = 496, ///< Location of the last state record
                    continuation = 504; ///< Location of the header continuation block
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
            static constexpr unsigned int size = HEADER::size; ///< Size of the continuation block
            static constexpr unsigned int states = (size / 8) - 1; ///< Number of states stored in a continuation block
            static constexpr int64_t next_cblock = size - 8; ///< Relative position of the pointer to the next continuation block
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

        /** Reads a state starting at the current file position.  The state data structure is as
         * follows:
         *
         *     u64      t (simulation time period)
         *     u32      numReaders
         *     u32      numBooks
         *     READER*numReaders    reader data
         *     BOOK*numBooks        book data
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
         *     u64              id
         *     dbl*DIM          position (DIM = dimensions)
         *     u64[]            friend ids
         *     (u8,u64,dbl)[]   library; see below for the u8 values.
         *     dbl              u
         *     dbl              u_lifetime
         *     dbl              cost_fixed
         *     dbl              cost_unit
         *     dbl              income
         *     BELIEF           profit belief
         *     BELIEF           profit extrapolated belief
         *     BELIEF           demand belief
         *     BELIEF           quality belief
         *     BELIEF[]         profit stream beliefs (for different K() values)
         *
         * where type[] indicates an array structured as:
         *     u32          length
         *     type*length  type
         *
         * and (type1,type2) indicates a single type1 value followed immediately by a single type2
         * value.
         *
         * The u8 in the library indicates the property of the book made up of the following bits:
         *     1 - set if this reader wrote the book
         *     2 - set if the reader pirated the book (instead of purchasing)
         *     4 - set if the book is new (i.e. added to library this period)
         *     >4 - reserved
         * 1 is exclusive of the other bits (they may not be set when 1 is set).
         *
         * BELIEF is a i64 belief location (or special value) and possibly a set of belief data, as
         * handled by readBelief(int64_t).
         *
         * profit stream beliefs may not be placeholder beliefs (i.e. default constructed objects);
         * such objects should simply be omitted when writing the data.  Each profit stream belief K
         * value must also be unique.
         */
        std::pair<eris::eris_id_t, ReaderState> readReader() const;

        /// Writes a reader state at the current file position.
        void writeReader(const ReaderState& reader);

        /// Structure holding parsed belief data
        typedef struct {
            uint32_t K; ///< Number of parameters; K=0 for a default constructed (invalid) model (the remaining values will be uninitialized)
            bool noninformative; ///< True if this is a noninformative model (in which case the remaining values will be uninitialized)
            Eigen::VectorXd beta; ///< belief::Linear beta vector
            double s2; ///< belief::Linear s2 value
            double n; ///< belief::Linear n value
            Eigen::MatrixXd V; ///< belief::Linear V matrix
        } belief_data;

        /** Reads a belief structure from the current file location.  The first field read is the
         * location value which has the following interpretations:
         *     - 0 indicates a default-constructed model (i.e. a placeholder object, not a real model)
         *     - a negative value in [-100,-1], indicates a noninformative model with `-val`
         *       parameters (e.g. -6 indicates a noninformative, 6-parameter model)
         *     - the value -512 indicates that the belief record immediately follows the current
         *       position
         *     - otherwise the value indicates the file position containing the belief record
         *       (identical belief records may be multiply referenced to save file space).  If the
         *       location is invalid (negative below -100, inside the file header, or past the end
         *       of the file) will raise an exception.
         *
         * In the case of a default-constructed model or a noninformative model, no further reading
         * is required.  When a file location (or the "immediately following" -512 value) is
         * recognized, that indicated location contains a belief record structured as follows:
         *
         *     u32      K, the number of model parameters
         *     dbl*K    beta vector (K values)
         *     dbl      s2
         *     dbl      n
         *     dbl*Z    the `Z=K(K+1)/2` elements of the lower triangle of the V matrix, in column-major
         *              order (that is, [0,0], [1,0], [2,0], ..., [K,0], [1,1], [2,1], ..., [K,1], ..., [K,K])
         *
         * The file's current position when this method returns will be the position of the next field.
         * the may be moved by this method and is not restored; if you
         * require restoration of the current position, handle it in the code calling this method.
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
        void writeBelief(const belief::Linear &belief);

        /** Structure for storing locations of previous beliefs.  The outer key stores K, the number
         *
         * of parameters of the model (which also precisely defines the record size).  The inner key
         * stores the (overflowed) sum of the belief record (not including the initial K parameter)
         * interpreted as an array of u32s, while the inner values are the file locations where the
         * record is stored.  Duplicates are possible: one must actually read the matched records to
         * verify whether they are exactly identical and not use ones that are different.
         *
         * This object only stores locations of beliefs that have been written or parsed; records
         * from a pre-existing file are not read (and thus unknown) until specifically accessed.
         */
        mutable std::unordered_map<uint32_t, std::unordered_multimap<uint32_t, int64_t>> belief_locations_;

        /** Returns the belief size (in terms of number of doubles) for an belief record of K
         * parameters, not including the initial u32 K value.  Multiple by `sizeof(double)` to get
         * the size in bytes. */
        static constexpr uint32_t beliefRecordSize(unsigned int K) {
            return
                K // beta
                + 2 // s2 and n
                + K*(K+1)/2; // V lower triangle
        }

        /** Reads a book state data from the current file position and returns it in a <eris_id_t,
         * BookState> pair, where the eris_id_t is the book id.  The book data is:
         *
         *     u64          id
         *     u64          author id
         *     dbl*DIM      position (DIM = dimensions)
         *     dbl          quality
         *     dbl          price (market is derived from this: market=true unless price is NaN)
         *     dbl          revenue
         *     dbl          revenueLifetime
         *     u64          sales
         *     u64          salesLifetime
         *     u64          pirated
         *     u64          piratedLifetime
         *     u64          copies
         *     u64          age
         *     u64          created
         *     u64          lifetime
         */
        std::pair<eris::eris_id_t, BookState> readBook() const;

        /// Whether the file contains a (non-default) state_begins_ value
        bool wrote_sharing_begins_ = false;

        /** Writes a book at the current file position.  See readBook() for data layout. */
        void writeBook(const BookState &book);

};

#if BOOST_ENDIAN_BIG_BYTE
#include "creativity/state/FileStorage-BE.hpp"
#elif BOOST_ENDIAN_LITTLE_BYTE
#include "creativity/state/FileStorage-LE.hpp"
#else
#error System endianness not supported (neither big-byte nor little-byte endianness detected)!
#endif

}}
