#include "creativity/data/CSVParser.hpp"
#include <sstream>

namespace creativity { namespace data {

CSVParser::CSVParser(const std::string &filename) {

    f_.exceptions(f_.failbit | f_.badbit);
    f_.open(filename, std::ios_base::in);

    std::string header;
    lineno_ = 1;
    std::getline(f_, header);
    header_ = split(header);

    updateFields();
}

bool CSVParser::readRow() {
    // Make sure there's at least one character available; if not, eof will be triggered.
    f_.peek();

    if (f_.eof()) return false;

    std::string row;
    std::getline(f_, row); // throws on error

    lineno_++;

    // Skip commented rows:
    if (row[0] == '#') return readRow();

    auto fields = split(row);
    if (fields.size() > header_.size() or header_.size() - fields.size() > allow_missing_values)
        throw std::invalid_argument("Invalid data on line " + std::to_string(lineno_) + ": number of fields (" + std::to_string(fields.size())
                + ") differs from that of the header (" + std::to_string(header_.size()) + ")");

    row_.clear();
    if (row_.capacity() < fields.size()) row_.reserve(fields.size());
    row_skipped_.clear();
    for (size_t fieldnum = 0; fieldnum < fields.size(); fieldnum++) {
        if (skip_.count(header_[fieldnum])) {
            row_skipped_[header_[fieldnum]] = fields[fieldnum];
        }
        else {
            auto &field = fields[fieldnum];
            size_t parse_pos = 0;
            double d;
            try {
                d = std::stod(field, &parse_pos);
                if (parse_pos != field.length()) throw std::invalid_argument("trailing garbage");
            }
            catch (const std::invalid_argument &e) {
                throw std::invalid_argument("CSVParser::readRow: found non-double value `" + field + "' for " +
                        header_[fieldnum] + ", line " + std::to_string(lineno_) + " (" + e.what() + ")");
            }
            row_.emplace_back(std::move(d));
        }
    }

    return true;
}

bool CSVParser::eof() const {
    return f_.eof();
}

CSVParser::iterator CSVParser::begin() {
    return iterator(*this, false);
}

CSVParser::iterator CSVParser::end() {
    return iterator(*this, true);
}

std::vector<std::string> CSVParser::split(const std::string &csr) {
    std::vector<std::string> results;
    std::stringstream ss(csr);
    std::string val;
    while (std::getline(ss, val, ',')) {
        results.push_back(val);
    }
    return results;
}

const std::unordered_set<std::string>& CSVParser::skip() const {
    return skip_;
}

void CSVParser::skip(const std::string &name) {
    skip_.insert(name);
    updateFields();
}

void CSVParser::dontSkip(const std::string &name) {
    skip_.erase(name);
    updateFields();
}

const std::vector<std::string>& CSVParser::header() const {
    return header_;
}

const std::vector<std::string>& CSVParser::fields() const {
    return fields_;
}

bool CSVParser::hasField(const std::string &field) const {
    if (field_pos_.count(field) > 0) return true;
    if (field_pos_.empty() and not fields_.empty()) {
        field_pos_.reserve(fields_.size());
        unsigned i = 0;
        for (auto f : fields_) {
            field_pos_.emplace(std::move(f), i++);
        }
        return field_pos_.count(field) > 0;
    }
    return false;
}

size_t CSVParser::fieldPosition(const std::string &field) const {
    if (field_pos_.empty() and not fields_.empty()) {
        hasField(field); // will populate field_pos_
    }
    return field_pos_.at(field);
}

double CSVParser::field(const std::string &field) const {
    if (row_.size() == 0) throw std::logic_error("Can't call field without a current row");

    return row_[fieldPosition(field)];
}

void CSVParser::updateFields() {
    fields_.clear();
    fields_.reserve(header_.size() - skip_.size());
    for (const auto &f: header_) {
        if (not skip_.count(f)) fields_.push_back(f);
    }
}

}}
