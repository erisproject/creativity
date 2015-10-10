#include "creativity/data/CSVParser.hpp"
#include "creativity/data/SUR.hpp"
#include "creativity/data/tabulate.hpp"
#include "creativity/cmdargs/Results.hpp"
#include <Eigen/Core>
#include <cerrno>
#include <exception>
#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <list>
#include <unordered_map>

using creativity::cmdargs::Results;
using namespace creativity::data;
using namespace eris;
using namespace Eigen;

// rawdata holds the data for rows with positive number of books written under piracy/public,
// raw_no_writing holds the data for rows with no books written under one or more situations.
MatrixXd rawdata, raw_no_writing;
std::unordered_map<std::string, std::shared_ptr<const SimpleVariable>> data_;
std::unordered_map<std::string, std::shared_ptr<const SimpleVariable>> data_nw_;
// Map row numbers to filenames:
std::vector<std::string> data_source, data_nw_source;
// Track fields that don't exist in some of the stages:
std::unordered_map<std::string, int> pre_nan_, piracy_nan_, public_nan_;
bool piracy_data = false, // True if we have piracy data
     public_data = false, // True if we have public data
     piracy_sr = false, // True if we have short-run piracy data
     public_sr = false; // True if we have short-run public sharing data

// Stupid wrapper class so that data["foo"] works (so data.at("foo") isn't needed)
class DataWrapper {
    public:
        DataWrapper(decltype(data_) &data) : d_{data} {}
        std::shared_ptr<const SimpleVariable> operator[](const std::string &field) {
            try { return d_.at(field); }
            catch (const std::out_of_range &re) {
                std::cerr << "Data exception: data field '" << field << "' does not exist\n";
                throw;
            }
        }
        const std::shared_ptr<const SimpleVariable>& operator[](const std::string &field) const {
            try { return d_.at(field); }
            catch (const std::out_of_range &re) {
                std::cerr << "Data exception: data field '" << field << "' does not exist\n";
                throw;
            }
        }
        bool has(const std::string &field) const { return d_.count(field) > 0; }
    private:
        decltype(data_) &d_;
};
DataWrapper data(data_), data_nw(data_nw_);


// Generate a row of data from a CSV row
void generateRow(
        const CSVParser &csv, // the CSV parser (current row will be used)
        Ref<RowVectorXd, 0, InnerStride<>> newrow, // Row to set
        const std::unordered_map<std::string, int> &data_column, // field to column map
        const std::unordered_map<std::string, int> &nans, // Fields to set to nan
        bool &no_writing, // Will be set true if books_written is less than 0.2
        const std::string &prefix, // Prefix to match (in addition to param.)
        const std::string &notprefix = "" // If non-empty, prefix to avoid even if prefix matches--typically `whatever.SR.`
        ) {

    for (unsigned fc = 0; fc < csv.fields().size(); fc++) {
        auto &d = csv.fields()[fc];
        if (
                d.substr(0, 6) == "param." or // Always include parameters
                (d.substr(0, prefix.size()) == prefix // Include things with the requested prefix
                    // ... unless `notprefix` is given and it matches the field
                    and not (!notprefix.empty() and d.substr(0, notprefix.size()) == notprefix))
           ) {
            newrow[data_column.at(d)] = csv.row()[fc];
        }
    }
    for (auto &nan : nans) newrow[nan.second] = std::numeric_limits<double>::quiet_NaN();
    if (piracy_data) newrow[data_column.at("piracy")] = prefix.substr(0, 7) == "piracy." ? 1 : 0;
    if (public_data) newrow[data_column.at("public")] = prefix.substr(0, 7) == "public." ? 1 : 0;
    if (piracy_sr or public_sr) newrow[data_column.at("shortrun")] = prefix.rfind(".SR.") == std::string::npos ? 0 : 1;

    if (
            (data_column.count(".books_written") and newrow[data_column.at(".books_written")] < 0.2)
        or
            (data_column.count(".book_p0") and std::isnan(newrow[data_column.at(".book_p0")]))
       )
        no_writing = true;
}



void readCSV(const std::string &filename) {
    CSVParser csv(filename);

    // The data contains values like pre.whatever, piracy.whatever, piracy.SR.whatever,
    // public.whatever, public.SR.whatever.  We need to convert those into five rows:
    // - one with `whatever' set to pre.whatever and piracy=public=shortrun=0
    // - one with `whatever' set to piracy.SR.whatever, piracy=1, public=0, shortrun=1
    // - one with `whatever' set to piracy.whatever, piracy=1, public=0
    // - one with `whatever' set to public.SR.whatever, piracy=0, public=1, shortrun=1
    // - one with `whatever' set to public.whatever, piracy=0, public=1
    //
    // The data might not, however, have any piracy and/or public and/or short run data, in which
    // case we omit the relevant row(s).

    // This maps every CSV field (by name) into a column.  The mapping isn't unique, however: any
    // pre.whatever, piracy.whatever, and public.whatever should all map to the same column, albeit
    // with different rows, as described above.
    //
    // For the various pre/piracy/public fields, we also store the common location in `_whatever`,
    // so that it's easier to look up.  Note: elements with leading single underscores should be
    // considered special!
    std::unordered_map<std::string, int> data_column;

    int next_col = 0;
    for (auto &f : csv.fields()) {
        if (f[0] == '.') throw std::runtime_error("CSV file has invalid field `" + f + "': fields may not start with .");
        std::string data_field;
        if (f.substr(0, 4) == "pre.") data_field = f.substr(3);
        else if (f.substr(0, 10) == "piracy.SR.") { data_field = f.substr(9); piracy_sr = true; }
        else if (f.substr(0, 10) == "public.SR.") { data_field = f.substr(9); public_sr = true; }
        else if (f.substr(0, 7) == "piracy.") { data_field = f.substr(6); piracy_data = true; }
        else if (f.substr(0, 7) == "public.") { data_field = f.substr(6); public_data = true; }

        if (not data_field.empty()) {
            if (data_column.count(data_field)) {
                // One of the pre/piracy/public equivalent values for this value was already seen, reuse that column
                data_column.insert({f, data_column[data_field]});
            }
            else {
                // This field hasn't been seen before
                data_column.insert({f, next_col});
                data_column.insert({data_field, next_col});
                next_col++;
            }
        }
        else {
            // Not a pre/piracy/public field, so just copy it as-is
            data_column.insert({f, next_col});
            next_col++;
        }
    }

    if (piracy_data) data_column.insert({"piracy", next_col++}); // piracy dummy
    if (public_data) data_column.insert({"public", next_col++}); // public sharing dummy
    if (piracy_sr or public_sr) data_column.insert({"shortrun", next_col++}); // shortrun periods dummy

    // Increase the matrix by 2520 row increments: 2520 is the lowest integer divisible by all
    // numbers from 1 to 10, so this will work properly for adding anywhere from 1 to 10 rows of
    // data per input row.  (Currently we have 1-5 possible rows, and 2520 is a reasonably sized
    // number for this, as well).
    constexpr unsigned rowincr = 2*2*2*3*3*5*7;
    rawdata.resize(0, next_col);
    raw_no_writing.resize(0, next_col);

    // Some fields don't have all three (for example, books_pirated has piracy_ and public_ values
    // but not a pre_ value; books_public_copies is only under public_), so we need to fill in some
    // NaNs for the ones where they don't exist: track the columns needing NaNs here:
    for (auto &pair : data_column) {
        if (pair.first[0] == '.') {
            if (not data_column.count("pre" + pair.first))
                pre_nan_.emplace(pair.first.substr(1), pair.second);
            if (piracy_data and not data_column.count("piracy" + pair.first))
                piracy_nan_.emplace(pair.first.substr(1), pair.second);
            if (public_data and not data_column.count("public" + pair.first))
                public_nan_.emplace(pair.first.substr(1), pair.second);
        }
    }

    unsigned raw_rows = 0, rawnp_rows = 0;
    while (csv.readRow()) {
        MatrixXd newrows(1 + piracy_data + public_data + piracy_sr + public_sr, rawdata.cols());
        unsigned rownum = 0;
        // Will be set to true if we find one or more periods with no writing/private marketing
        bool no_writing = false;

        // First row: pre values
        generateRow(csv, newrows.row(rownum++), data_column, pre_nan_, no_writing, "pre.");

        // Second: short run piracy
        if (piracy_sr) {
            generateRow(csv, newrows.row(rownum++), data_column, piracy_nan_, no_writing, "piracy.SR.");
        }
        // Third: long run piracy
        if (piracy_data) {
            generateRow(csv, newrows.row(rownum++), data_column, piracy_nan_, no_writing, "piracy.", "piracy.SR.");
        }
        // Fourth: short run public
        if (public_sr) {
            generateRow(csv, newrows.row(rownum++), data_column, public_nan_, no_writing, "public.SR.");
        }
        // Fifth: long run public
        if (public_data) {
            generateRow(csv, newrows.row(rownum++), data_column, public_nan_, no_writing, "public.", "public.SR.");
        }

        MatrixXd &rawmatrix = (no_writing ? raw_no_writing : rawdata);
        auto &rawrows = (no_writing ? rawnp_rows : raw_rows);
        auto &source = (no_writing ? data_nw_source : data_source);
        if (rawrows + newrows.rows() >= rawmatrix.rows()) {
            rawmatrix.conservativeResize(rawmatrix.rows() + rowincr, NoChange);
        }
        rawmatrix.middleRows(rawrows, newrows.rows()) = newrows;
        unsigned need = rawrows + newrows.rows();
        if (source.size() < need) source.resize(need);
        for (unsigned i = 0; i < newrows.rows(); i++) source[rawrows + i] = csv.rowSkipped().at("source");
        rawrows += newrows.rows();
    }

    // Resize the probably-oversized data matrices back to the number of rows we actually filled
    if (rawdata.rows() > raw_rows) rawdata.conservativeResize(raw_rows, NoChange);
    if (raw_no_writing.rows() > rawnp_rows) raw_no_writing.conservativeResize(rawnp_rows, NoChange);

    for (auto &dc : data_column) {
        if (dc.first.substr(0, 4) == "pre." or dc.first.substr(0, 7) == "public." or dc.first.substr(0, 7) == "piracy.")
            continue;
        std::string name(dc.first);
        if (dc.first[0] == '.') name = dc.first.substr(1);
        data_.insert({name, SimpleVariable::create(name, rawdata.col(dc.second))});
        data_nw_.insert({name, SimpleVariable::create(name, raw_no_writing.col(dc.second))});
    }
}

int main(int argc, char *argv[]) {
    Results args;
    args.parse(argc, argv);

    try {
        readCSV(args.input);
    }
    catch (const std::exception &e) {
        std::cerr << "Error: failed to read input file `" << args.input << "': " << e.what() << "\n\n";
    }

    unsigned nobs_per_sim = 1 + piracy_data + public_data + public_sr + piracy_sr;
    unsigned nobs_w_writing = data["net_u"]->size(), nobs_wo_writing = data_nw["net_u"]->size();
    unsigned sims_w_writing = nobs_w_writing / nobs_per_sim, sims_wo_writing = nobs_wo_writing / nobs_per_sim;
    std::cout << "Data summary:\n" <<
        "    " << nobs_w_writing + nobs_wo_writing << " total observations (from " << sims_w_writing+sims_wo_writing << " simulations)\n" <<
        "    " << nobs_w_writing << " observations (from " << sims_w_writing << " simulations) with non-zero # books written during each situation\n" <<
        "    " << nobs_wo_writing << " observations (from " << sims_wo_writing << " simulations) with zero books written during one or more situations\n";


    SUR avg_effects;
    for (auto &y : {"net_u", "books_written", "book_quality", "book_p0", "book_revenue", "book_profit"}) {
        Equation eq(data[y]);
        eq % 1;
        if (piracy_data) eq % data["piracy"];
        if (public_data) eq % data["public"];
        avg_effects.add(std::move(eq));
    }
    avg_effects.gather();
    for (unsigned i = 0; i < avg_effects.equations().size(); i++) {
        auto yi = avg_effects.y(i);
        for (unsigned r = 0; r < yi.size(); r++) {
            if (std::isnan(yi[r])) {
                std::cerr << "Error: found NaN for " << avg_effects.equations()[i].depVar()->name() << ", source file " << data_source[r] << "\n";
            }
        }
    }

    //std::cout << tabulate(avg_effects);
    
    avg_effects.solve();
    std::cout << "Average effects:\n================\n" << avg_effects;

    SUR marg_effects;
    for (auto &y : {"net_u", "books_written", "book_quality", "book_p0", "book_revenue", "book_profit"}) {
        Equation eq(data[y]);
        eq % 1;
        if (piracy_data) eq % data["piracy"];
        if (public_data) eq % data["public"];
        for (auto &x : {"param.density", "param.cost_market", "param.cost_unit", "param.creation_time", "param.creation_fixed"}) {
            if (not data.has(x)) continue;
            if (not pre_nan_.count(x)) eq % data[x];
            if (piracy_data and not piracy_nan_.count(x)) eq % (data["piracy"] * data[x]);
            if (public_data and not public_nan_.count(x)) eq % (data["public"] * data[x]);
        }
        // These don't get included in "pre":
        for (auto &x : {"param.cost_piracy", "param.piracy_link_proportion"}) {
            if (not data.has(x)) continue;
            if (piracy_data and not piracy_nan_.count(x)) eq % (data["piracy"] * data[x]);
            if (public_data and not public_nan_.count(x)) eq % (data["public"] * data[x]);
        }
        // These don't get included in "pre" or "piracy":
        for (auto &x : {"param.public_sharing_tax"}) {
            if (not data.has(x)) continue;
            if (public_data and not public_nan_.count(x)) eq % (data["public"] * data[x]);
        }
        marg_effects.add(eq);
    }

    marg_effects.solve();
    std::cout << "\n\n\nMarginal effects:\n=================\n" << marg_effects;
}

