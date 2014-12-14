#pragma once
#include "creativity/state/StorageBackend.hpp"
#include <vector>
// Hack to make libpqxx not try to load <tr1/memory>, which isn't needed or wanted
#include <pqxx/compiler-public.hxx>
#undef PQXX_TR1_HEADERS
#undef PQXXTR1
#define PQXXTR1 std
// End hack
#include <pqxx/pqxx>
#include <iomanip>

namespace creativity { namespace state {

/** Class for PostgreSQL-based storage.  Note that most of the methods of this class will throw
 * errors when the underlying libpqxx library encounters an error.  You should not attempt to use
 * the PsqlStorage object after such an exception is thrown as the result is unpredictable: the
 * object could fail to write new data or could write data inconsistently.
 *
 * This object is also not thread-safe on its own: you must ensure that only one thread at a time
 * accesses it, even for read-only methods.
 *
 * See database.pgsql included in the creativity distribution for the required table structure.
 */
class PsqlStorage : public StorageBackend {
    public:
        PsqlStorage() = delete;

        /** Constructs and returns a PsqlStorage object that connects to the given database for
         * reading and writing state data.
         *
         * \param connect the libpqxx connection string specifying the database.  Anything supported
         * by libpq should work here (see
         * http://www.postgresql.org/docs/9.3/static/libpq-connect.html#LIBPQ-CONNSTRING)
         *
         * Some example connect strings:
         *
         *     "host=localhost port=5432 dbname=mydb connect_timeout=10"
         *     "host=somewhere username=myuser dbname=mydb password=mypassword requiressl=true"
         *     "postgresql://other@localhost:5432/otherdb?connect_timeout=10&application_name=myapp"
         *
         * \param id if passed as a non-zero integer, the given simulation will be read from the
         * database, and an std::out_of_range exception thrown if no simulation with that ID exists.
         * If the record was found, `.seed` is updated with the stored seed value.
         *
         * Otherwise, with the default value of 0, a new simulation record will be created and
         * become available in 'id' during construction.  The current value of eris::Random::seed()
         * is stored in the new records seed field, and also stored in `.seed`.
         *
         * \throws pqxx::pqxx_exception if a database error occurs
         * \throws std::out_of_range if the given id doesn't exist
         * \throws std::logic_error if the given id is invalid (i.e. negative but not -1)
         */
        PsqlStorage(const std::string &connect, unsigned int id = 0);

        /** Loads the requested state data from the database into a State object and returns it
         * (wrapped in a shared_ptr).
         */
        virtual std::shared_ptr<const State> load(eris::eris_time_t t) const override;

        /** Stores the given settings in the database.  If existing settings are present, they are
         * replaced.
         */
        void writeSettings(const CreativitySettings &settings) override;

        /** Reads the settings from the database and sets them into the given settings object.
         */
        void readSettings(CreativitySettings &settings) const override;

        /// Returns the number of states currently stored in the database.
        virtual size_t size() const override;

        /// The simulation id, as given or established during construction.
        const int32_t &id = id_;

        /// The simulation seed, as read or determined during construction.
        const int64_t &seed = seed_;

        /// Accesses the pqxx::connection object, locking it out from other threads.  Be careful.
        std::pair<pqxx::connection&, std::unique_lock<std::mutex>> connection();

    protected:
        /// Inserts the given state into the database.
        virtual void thread_insert(std::shared_ptr<const State> &&s) override;

    private:
        /// The simulation id (mutable but protected)
        int32_t id_;

        /// The simulation seed (mutable but protected);
        int64_t seed_;

        /// The number of dimensions (as last read or written to settings)
        mutable uint32_t dimensions_;

        /** The file buffer object. Mutable because we need to read from it in const methods. */
        mutable std::unique_ptr<pqxx::connection> conn_;

        /** Mutex guarding conn_. */
        mutable std::mutex conn_mutex_;

        /// Sets up various prepared queries, sets required settings, etc.
        void initialConnection();

        /// Reads a reader, returns it in a pair suitable for moving into State.readers.
        std::pair<eris::eris_id_t, ReaderState> readReader(const pqxx::tuple &reader_row, pqxx::work &trans) const;

        /// Adds a reader to the database associated with the given state id.
        void insertReader(const ReaderState &reader, int32_t state, pqxx::work &trans);

        /// Inserts a belief record for a reader
        void insertBelief(eris::eris_id_t dbid, const std::string &type, const belief::Linear &belief, pqxx::work &trans);

        /// Reads a book, returns it in a pair suitable for moving into State.books.
        std::pair<eris::eris_id_t, BookState> readBook(const pqxx::tuple &book_row, pqxx::work &trans) const;

        /** Parses a PG string of doubles such as `"{12.3,-4.567e65,NaN,-Infinity}"` into a vector
         * of doubles.
         *
         * \param doubles the string containing the PostgreSQL array of numeric values (typically
         * double precisions)
         * \param expected_size if provided, specifies the expected length of the array and raises
         * an exception if the actual size does not match.  Defaults to -1, which does not check.
         *
         * \throws std::invalid_argument for various parse conditions: if the string doesn't begin
         * and end with { and }; if any of the arguments could not be parsed as doubles; or if the
         * number of doubles doesn't match `expected_size`.
         *
         * \returns a vector of parsed double values.
         */
        static std::vector<double> parseDoubleArray(const std::string &doubles, int expected_size = -1);

        /** Converts a sequence of numeric values into a PostgreSQL array string of double precision
         * values.  Essentially this is the inverse of parseDoubleArray().
         *
         * \param start an iterator to the beginning of the sequence; typically container.begin()
         * \param end an iterator to the end of the squence; typically container.end()
         * \returns the string in the format expected by PostgreSQL
         */
        template <class Iterator>
        static typename std::enable_if<std::is_arithmetic<typename std::iterator_traits<Iterator>::value_type>::value, std::string>::type
        createDoubleArray(Iterator start, const Iterator &end) {
            std::ostringstream ss;
            ss << "{";
            // Set the significant digits to the max we need to represent the given type.  The
            // following gives 0 for integer types, and the maximum we need for float/double/long
            // double (which is 9 for float, 17 for double, and 36 for long double).  For things
            // giving 0, the setprecision is irrelevant anyway.
            auto precision = std::numeric_limits<typename std::iterator_traits<Iterator>::value_type>::max_digits10;
            if (precision > 0) ss.precision(precision);

            bool first = true;
            while (start != end) {
                if (first) first = false;
                else ss << ",";
                ss << *start;
                start++;
            }
            ss << "}";
            return ss.str();
        }
};

}}
