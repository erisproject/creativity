#include <creativity/data/Treatment.hpp>

namespace creativity { namespace data {

using namespace Eigen;

Treatment::Treatment(CSVParser &&csv) {
    readCSV(std::move(csv));
}

Treatment::Treatment(const std::string &filename) : Treatment(CSVParser(filename)) {}

void Treatment::readCSV(CSVParser &&csv) {
    // The data contains values like pre.whatever, piracy.whatever, piracy.SR.whatever,
    // public.whatever, public.SR.whatever.  We need to convert those into five rows:
    // - one with `whatever' set to pre.whatever and piracy=public=SR=0, LR=1
    // - one with `whatever' set to piracy.SR.whatever, piracy=1, public=0, SR=1, LR=0
    // - one with `whatever' set to piracy.whatever, piracy=1, public=0, SR=0, LR=1
    // - one with `whatever' set to public.SR.whatever, piracy=0, public=1, SR=1, LR=0
    // - one with `whatever' set to public.whatever, piracy=0, public=1, SR=0, LR=1
    //
    // The data might not, however, have any piracy and/or public and/or short run data, in which
    // case we omit the relevant row(s).
    //
    // This maps every CSV field (by name) into a column.  The mapping isn't unique, however: any
    // pre.whatever, piracy.whatever, and public.whatever all map to "whatever", albeit with
    // different rows, as described above.
    //
    // Fields starting with "param." are left as-is; anything else must have one of the "prefix."
    // values described above; anything else is unrecognized and will throw an exception.  Only the
    // "param." prefixes are preserved: all other prefixes are removed.

    if (have_data_) throw std::logic_error("Treatment::readCSV() can only be called once");
    have_data_ = true;

    csv.skip("source");

    int next_col = 0;
    for (auto &f : csv.fields()) {
        std::string data_field;
        if (f.substr(0, 4) == "pre.") { data_field = f.substr(4); has_pre_ = true; }
        else if (f.substr(0, 10) == "piracy.SR.") { data_field = f.substr(10); has_piracy_sr_ = true; }
        else if (f.substr(0, 10) == "public.SR.") { data_field = f.substr(10); has_public_sr_ = true; }
        else if (f.substr(0, 7) == "piracy.") { data_field = f.substr(7); has_piracy_ = true; }
        else if (f.substr(0, 7) == "public.") { data_field = f.substr(7); has_public_ = true; }
        else if (f.substr(0, 6) == "param.") { data_field = f; } // leave leading param.
        else { throw std::runtime_error("CSV file has invalid/unknown field `" + f + "': fields must have a (known) prefix"); }

        if (data_field.substr(0, 6) != "param." and data_field.find('.') != data_field.npos) {
            throw std::runtime_error("CSV file has invalid field `" + f + "': non-parameter fields may not contain . (except as field type prefixes)");
        }
        if (data_column_.count(data_field) == 0) {
            // It's a field that we haven't seen before: add a new column
            data_column_.emplace(data_field, next_col);
            next_col++;
        }
    }

    // The pre-data plus each treatment effect contributes an output row
    rows_per_sim_ = has_pre_ + has_piracy_ + has_piracy_sr_ + has_public_ + has_public_sr_;

    if (rows_per_sim_ == 0) throw std::runtime_error("CSV file contains no usable data (no pre, piracy, or public observations found)");

    // Check required data
    requirePre(require_pre_);
    requirePiracy(require_piracy_);
    requirePublic(require_public_);
    requireSR(require_sr_);

    // Add dummy columns (even if the data means they will always be 0/1)
    data_column_.insert({"piracy", next_col++}); // piracy dummy
    data_column_.insert({"public", next_col++}); // public sharing dummy
    data_column_.insert({"SR", next_col++}); // short-run stage dummy
    data_column_.insert({"LR", next_col++}); // long-run stage dummy

    // Increment the matrix by 1000 *file* rows at a time (i.e. we need to do a relatively expensive
    // matrix resize every 1000 input rows).
    const unsigned rowincr = 1000 * rows_per_sim_;
    data_ = MatrixXd(rowincr, next_col);

    unsigned pos = 0;
    while (csv.readRow()) {
        if (pos >= data_.rows())
            data_.conservativeResize(data_.rows() + rowincr, NoChange);

        // NB: this order matters!  TreatmentFilter relies on it.
        // First row: pre values
        generateRow(csv, data_.row(pos++), "pre.");
        // Second: short run piracy
        if (has_piracy_sr_) generateRow(csv, data_.row(pos++), "piracy.SR.");
        // Third: long run piracy
        if (has_piracy_) generateRow(csv, data_.row(pos++), "piracy.");
        // Fourth: short run public
        if (has_public_sr_) generateRow(csv, data_.row(pos++), "public.SR.");
        // Fifth: long run public
        if (has_public_) generateRow(csv, data_.row(pos++), "public.");

        source_.push_back(csv.rowSkipped().at("source"));
    }

    // Resize the probably-oversized data matrix back to the number of rows we actually filled
    if (data_.rows() > pos) data_.conservativeResize(pos, NoChange);
}

const std::unordered_map<std::string, unsigned>& Treatment::columns() const {
    return data_column_;
}

std::shared_ptr<const SimpleVariable> Treatment::variable(const std::string &name) const {
    if (var_cache_.count(name)) return var_cache_[name];

    unsigned col = data_column_.at(name); // Throws if not found
    // Insert and return:
    return var_cache_.emplace(name, SimpleVariable::create(name, data_.col(col))).first->second;
}


// Generate a row of data from a CSV row
void Treatment::generateRow(
        const CSVParser &csv,
        Ref<RowVectorXd, 0, InnerStride<>> newrow,
        const std::string &prefix) const {

    for (const auto &field : data_column_) {
        // Look for special dummies:
        if (field.first == "piracy") newrow[field.second] = (prefix.substr(0, 7) == "piracy.");
        else if (field.first == "public") newrow[field.second] = (prefix.substr(0, 7) == "public.");
        else if (field.first == "LR") newrow[field.second] = (prefix.find(".SR.") == prefix.npos);
        else if (field.first == "SR") newrow[field.second] = (prefix.find(".SR.") != prefix.npos);
        else {
            std::string csv_field = field.first.substr(0, 6) == "param."
                ? field.first
                : prefix + field.first;
            newrow[field.second] =
                csv.hasField(csv_field) ? csv.field(csv_field) : std::numeric_limits<double>::quiet_NaN();
        }
    }
}

#define REQUIRE(Meth, var, name) \
void Treatment::require##Meth(bool req) { \
    require_##var##_ = req; \
    if (have_data_ and require_##var##_ and not has_##var##_) throw std::runtime_error("CSV file is missing required `" name "' data"); \
}
REQUIRE(Pre, pre, "pre-piracy")
REQUIRE(Piracy, piracy, "piracy")
REQUIRE(Public, public, "public")
void Treatment::requireSR(bool req) {
    require_sr_ = req;
    if (have_data_ and require_sr_) {
        if (has_piracy_ and not has_piracy_sr_) throw std::runtime_error("CSV file has long-run but no short-run piracy data");
        if (has_public_ and not has_public_sr_) throw std::runtime_error("CSV file has long-run but no short-run public data");
    }
}

}}
