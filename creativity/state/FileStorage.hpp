#pragma once
#include "creativity/state/StorageBackend.hpp"
#include "creativity/BookCopy.hpp"
#include "creativity/CreativitySettings.hpp"
#include <eris/belief/BayesianLinearRestricted.hpp>
#include <eris/types.hpp>
#include <eris/serialize/serializer.hpp>
#include <eris/serialize/Serialization.hpp>
#include <Eigen/Core>
#include <boost/detail/endian.hpp>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <vector>
#include <iostream>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <map>
#include <utility>

namespace creativity { namespace state {

/** Class for file-based storage.  Note that most of the methods of this class will throw errors
 * when the underlying file stream encounters an error.  You should not attempt to use the
 * FileStorage object after such an exception is thrown as the result is unpredictable: the object
 * could simply stop writing new data to the file or the file could become corrupted.
 *
 * Writing new states is done by a background thread which is spawned the first time enqueue() is
 * called with a new state.
 */
class FileStorage : public StorageBackend, public eris::serialize::Serialization {
    public:
        /** The supported record types within a state. */
        enum TYPE : uint8_t {
            TYPE_DONE = 0, ///< Psuedo-type indicating the end of the record list
            TYPE_READERS = 1, ///< An array of readers
            TYPE_BOOKS = 2, ///< An array of books
            TYPE_PUBLIC_TRACKER = 3 ///< A PublicTracker state
        };

        /** Constructs a FileStorage that stores file content in an in-memory buffer.  Beware: this
         * usage requires considerably more memory for running a simulation, and the results are not
         * automatically stored.
         *
         * \param settings the CreativitySettings reference which is read when creating a new file
         * and updated when reading an existing file.
         */
        FileStorage(CreativitySettings &settings) : StorageBackend(settings) {
            memory();
        }

        /** Constructs and returns a FileStorage object that uses the given file for reading and
         * (optionally) writing state data.  The file is read or created immediately.  All arguments
         * (except the first) are forwarded to Serialize::open().
         *
         * \sa Serialize::open()
         *
         * \throws various exceptions if the file does not exist, cannot be read, is empty, or
         * contains invalid data.
         */
        template <typename... Args>
        FileStorage(CreativitySettings &settings, const std::string &filename, Mode mode, Args&&... args)
        : StorageBackend(settings)
        {
            open(filename, mode, std::forward<Args>(args)...);
        }

        /** Default move constructor. */
        //FileStorage(FileStorage&&) = default;

        /** Default move assignment operator. */
        //FileStorage& operator=(FileStorage&&) = default;

        /** Throws a eris::serialize::Serializer::parse_error exception.  The given message is
         * prefixed with `Parsing file failed [pos=123]: `, where `123` is the current file
         * position.
         */
        void throwParseError(const std::string& message) const;

        /** Loads the requested state data from the open file into a State object and returns it
         * (wrapped in a shared_ptr).
         */
        virtual std::shared_ptr<const State> load(eris::eris_time_t t) override;

        /// Returns the number of states currently stored in the file.
        virtual size_t size() override;

        /** Flushes the file stream.  This calls `flush()` on the underlying file object.  This is
         * normally not required: any changes will be automatically flushed when the FileStorage
         * object is destroyed; calling this manually primarily guards against unflushed data being
         * unwritten if the process is terminated abnormally (i.e. via a signal or a segfault).
         *
         * \sa creativity::state::StorageBackend::flush
         * \sa std::basic_ostream::flush
         */
        virtual void storage_flush() override;

    protected:

        /// Registers the various creativity parameters as header fields.
        void configureHeaderFields() override;

        /// Calls updateHeaderFields() to rewrite settings
        void writeSettings() override;

        /// Adds a state pointer block to the end of the header.
        void writeExtraHeader() override;

        /// Stores the location of the state pointer block at the end of the head.
        void readExtraHeader() override;

    private:

        /// The location of the state pointer and main reader/libraries list blocks
        uint64_t state_pointer_block_, readerlib_block_;

        // Sorting method that sorts by the library item's acquired date
        class lib_comp_less final {
            public:
                using ref = const std::pair<const eris::eris_id_t, BookCopy>&;
                bool operator()(ref a, ref b) const {
                    return a.second.acquired < b.second.acquired;
                }
        };

        // Reader library data
        struct lib_data {
            // book id to BookCopy
            std::map<uint64_t, BookCopy> library;
            // library, but sorted by acquired date (for fast retrieval)
            std::multiset<std::reference_wrapper<std::pair<const uint64_t, BookCopy>>, lib_comp_less> library_acq_sorted;
        };

        /// Cached individual reader library data already stored in the file; the key is the reader
        //id, the value is the pair of / the block location and the library data content.
        std::map<eris::eris_id_t, std::pair<uint64_t, lib_data>> reader_lib_;

        /** The file version (1).  Note that there were earlier FileStorage versions that were not
         * based on eris::Serializer at all, and have an incompatible file format; this 1 refers to
         * the first Serializer-compatible version.
         */
        virtual uint32_t appFileVersion() const override { return 1; }

        /** The file application name, "creativity". */
        virtual std::string appName() const override { return "creativity"; }

        /// Called from the queue thread to write the given State to the file.
        virtual void thread_insert(std::shared_ptr<const State> &&s) override;

            /** a reader library record size: the book id (u32), the acquisition period (u32), the
             * reader-specific quality (dbl), and the library book status (u8: 0=wrote, 1=bought,
             * 2=pirated).
             */

        /** Ensures that the library reader list contains all readers, and that each reader's
         * library block contains all of the reader's library books, adding anything as needed.
         */
        void updateLibraries(const std::map<eris::eris_id_t, ReaderState> &r);

        /** Reads a simulation state starting at the current file position.  The state data
         * structure is as follows:
         *
         *     u32      t (simulation time period)
         *
         * followed by any number of:
         *     u8       TYPE != 0
         *     DATA     specific to TYPE
         *
         * and finally a terminating u8 value of TYPE = 0.  Typically there is at least a reader
         * array (TYPE_READERS) and (except in the first states) a book array (TYPE_BOOKS).
         *
         * \see TYPE for the different TYPE values supported.
         */
        std::shared_ptr<const State> readState();

        /// Writes the given simulation state at the current file position.
        void writeState(const State &state);

        /** Reads a ReaderState record from the current file position and returns it in an
         * {eris_id_t, ReaderState} pair, where `.first` is the id.  Such a record consists of:
         *
         *     u32              id
         *     dbl*DIM          position (DIM = dimensions)
         *     u32[]            friend ids
         *     dbl              u
         *     dbl              u_lifetime
         *     dbl              creation_shape
         *     dbl              creation_scale
         *     BELIEF           profit belief
         *     BELIEF           profit extrapolated belief
         *     BELIEF           demand belief
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
         * BELIEF is a set of belief data, as handled by readBelief().  profit extrapolated belief
         * is set to a belief only if it differs from profit_belief; otherwise it is simply set to a
         * no-data, noninformative belief record (and shouldn't be used).
         *
         * profit stream beliefs may not be placeholder beliefs (i.e. default constructed objects);
         * such objects should simply be omitted when writing the data.  Each profit stream belief K
         * value must also be unique.
         */
        std::pair<eris::eris_id_t, ReaderState> readReader(eris::eris_time_t t);

        /// Writes a reader state at the current file position.
        void writeReader(const ReaderState& reader);

        /// Structure holding parsed belief data
        typedef struct {
            uint32_t K = 0; ///< Number of parameters; K=0 for a default constructed (invalid) model (the remaining values will be uninitialized)
            bool noninformative = true; ///< True if this is a noninformative model (in which case the following are not set)
            Eigen::VectorXd beta; ///< eris::belief::BayesianLinear beta vector
            double s2; ///< eris::belief::BayesianLinear s2 value
            double n; ///< eris::belief::BayesianLinear n value
            Eigen::MatrixXd Vinv; ///< eris::belief::BayesianLinear Vinv matrix
            eris::belief::BayesianLinearRestricted::DrawMode last_draw_mode; ///< The last draw mode from this belief (for restricted models)
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
         * - bit 2 (bit & 2): if set on a restricted model, this indicates the last draw from this
         *   model used Gibbs sampling (if not a restricted model, this bit is unused).
         * - bit 3 (bit & 4): the model is noninformative.  If this bit is set, then instead of
         *   beta/V/etc. data the record contains non-informative data that hasn't been loaded yet;
         *   otherwise the beta/V/n/s2 data follows.
         * - bit 4 (bit & 8): if set on a restricted model, this indicates the last draw from this
         *   model used rejection sampling (if not a restricted model, this bit is unused).  If
         *   neither this bit nor bit 2 are set, the last draw mode was Auto (which means no draw
         *   took place).
         * - other bits are currently unused.
         */
        belief_data readBelief();

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
         *     u8           status fields (current just 1 for private market)
         *     dbl          price (market is derived from this: market=true unless price is NaN)
         *     dbl          revenue
         *     dbl          revenue_lifetime
         *     dbl          prize
         *     dbl          prize_lifetime
         *     u32          sales
         *     u32          sales_lifetime_private
         *     u32          sales_lifetime_public
         *     u32          pirated
         *     u32          piratedLifetime
         *     u32          created
         *     u32          lifetime_private
         */
        std::pair<eris::eris_id_t, BookState> readBook();

        /** Writes a book at the current file position.  See readBook() for data layout. */
        void writeBook(const BookState &book);

        /** Reads public tracker state data from the current file position and returns it in a
         * PublicTrackerState unique pointer.  The data is:
         *
         *     u32          id
         *     dbl          tax (per-user lump sum tax)
         *     dbl          unspent (!= 0 only if the previous period had no publid downloads)
         */
        std::unique_ptr<PublicTrackerState> readPublicTracker();

        /// Writes a public tracker to the current file position.
        void writePublicTracker(const PublicTrackerState &pt);
};

}}
