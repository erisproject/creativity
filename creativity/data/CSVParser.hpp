#pragma once
#include <eris/noncopyable.hpp>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <Eigen/Core>

namespace creativity { namespace data {

/** Primitive comma-separated-value file parser.  This is far from a complete CSV parser: it doesn't
 * handle non-numeric values at all, for example; it is simply intended to parse the numeric CSV
 * output produced by data.cpp.
 *
 * The intended basic usage is:
 *
 *     CSVParser csv("filename.csv");
 *     for (auto &row : csv) {
 *         // Do something with row (which is a `const Eigen::RowVectorXd`)
 *     }
 *
 *     or, equivalently:
 *
 *     CSVParser csv("filename.csv");
 *     while (csv.readRow()) {
 *         // Do something with csv.row()
 *     }
 *
 * One non-standard extension is that this allows comment lines beginning with # in the file: any
 * such line will be skipped.
 */
class CSVParser : private eris::noncopyable {
    public:
        /// Not default constructible
        CSVParser() = delete;
        /** Opens a csv file for parsing.
         *
         * \param filename the file to read
         * \throws std::ios_base::failure for underlying IO errors
         */
        explicit CSVParser(const std::string &filename);

        /** The set of fields to skip; all non-numeric fields must be added here before attempting
         * to parse a line.  The default contains 'source'.
         *
         * \sa skip(const std::string &name)
         * \sa dontSkip()
         */
        const std::unordered_set<std::string>& skip() const;

        /** Adds the given field name to the list of header fields to skip, if not already present. */
        void skip(const std::string &name);

        /** Removes the given field name from the list of header fields to skip, if present. */
        void dontSkip(const std::string &name);

        /** Returns the vector of header names read during construction.  This includes fields that
         * will be skipped; if you don't want them, call fields() instead.
         *
         * \sa fields()
         */
        const std::vector<std::string>& header() const;

        /** Returns the vector of field names corresponding to a readRow() call.  This is
         * essentially the same as header(), but skips any fields in skip().
         *
         * \sa header()
         * \sa skip()
         */
        const std::vector<std::string>& fields() const;

        /** Reads the next line of the CSV file, storing it in row(), replacing what was previously
         * stored there.  Returns true if a row was read, false if the end of the file was hit,
         * throws an exception if something goes wrong.
         *
         * \throws std::ios_base::failure if a read error occurs
         * \throws std::logic_error if readRow() is called when eof() is true (i.e. the previous
         * readRow() already hit the end of the file and returned false).
         * \throws std::invalid_argument if one of the fields to parse (i.e. all of those not in
         * skip()) cannot be converted to a double value.
         */
        bool readRow();

        /** Returns true if the parser has reached the end of the file. */
        bool eof() const;

        /** Accesses the most-recently-read row of the file.  If no row has yet been read, this will
         * be an empty vector.
         */
        const Eigen::RowVectorXd& row() const { return row_; }

        /** Accesses any skipped fields in the most-recently-read row. */
        const std::unordered_map<std::string, std::string>& rowSkipped() const { return row_skipped_; }

        // forward declaration
        class iterator;

        /// Returns an iterator that reads through the file
        iterator begin();

        /// Returns a past-the-end iterator
        iterator end();

    private:
        // Split a string by , and return a vector of elements
        static std::vector<std::string> split(const std::string &csr);

        void updateFields(); // regenerate fields, omitting things in skip_
        std::unordered_set<std::string> skip_{{"source"}}; // things to skip in header_ when making fields_
        std::vector<std::string> header_; // the header read from the file
        std::vector<std::string> fields_; // the header fields remaining after skipping skip_
        std::fstream f_;
        size_t lineno_; // Tracks the current line number
        Eigen::RowVectorXd row_; // The most-recently-read row (reused)
        std::unordered_map<std::string, std::string> row_skipped_; // Any skipped fields on most-recently-read row
};

/** Iterator class that allows iterating through the file.  The iterator satisfies the requirements
 * of a ForwardIterator. */
class CSVParser::iterator final : public std::iterator<std::input_iterator_tag, const Eigen::RowVectorXd, long> {
    public:
        /// Dereferences the iterator, returning the current CSVParser row.
        reference operator*() { return csv_.row(); }

        /// Dereferences the iterator, returning the current CSVParser row pointer.
        pointer operator->() { return &csv_.row(); }

        /** Increments the iterator, reading the next row of the file.  Previous iterators are
         * invalidated.
         */
        iterator& operator++() { end_ = not csv_.readRow(); return *this; }

        /** Return true if the given object is a reference to the current object, or if the current
         * object is at the end of file and the given object is a special end iterator for the same
         * CSVParser object.
         */
        bool operator==(const iterator &other) {
            return
                &csv_ == &(other.csv_) // same CSV object
                and
                end_ == other.end_; // Both are end (either by calling end() or by hitting eof)
        }

        /** Returns the negation of the == operator. */
        bool operator!=(const iterator &other) { return !(*this == other); }

    private:
        iterator() = delete;
        CSVParser &csv_;
        bool end_; // true if this is a past-the-end iterator
        iterator(CSVParser &csv, bool end) : csv_(csv), end_(end) { if (!end_) csv_.readRow(); }
        friend class CSVParser;
};

}}
