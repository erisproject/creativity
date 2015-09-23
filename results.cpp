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

Eigen::MatrixXd rawdata;
std::unordered_map<std::string, std::shared_ptr<const SimpleVariable>> data_;
// Track fields that don't exist in some of the stages:
std::unordered_map<std::string, int> pre_nan_, piracy_nan_, public_nan_;
bool piracy_data = false, public_data = false;
// Stupid wrapper class so that data["foo"] works (so data.at("foo") isn't needed)
class {
    public:
        std::shared_ptr<const SimpleVariable> operator[](const std::string &field) {
            try { return data_.at(field); }
            catch (const std::out_of_range &re) {
                std::cerr << "Data exception: data field '" << field << "' does not exist\n";
                throw;
            }
        }
        const std::shared_ptr<const SimpleVariable>& operator[](const std::string &field) const {
            try { return data_.at(field); }
            catch (const std::out_of_range &re) {
                std::cerr << "Data exception: data field '" << field << "' does not exist\n";
                throw;
            }
        }
        bool has(const std::string &field) const { return data_.count(field) > 0; }
} data;

void readCSV(const std::string &filename) {
    CSVParser csv(filename);

    // The data contains values like pre.whatever, piracy.whatever, public.whatever.  We need to convert those into three rows:
    // - one with whatever set to pre.whatever and piracy=public=0
    // - one with whatever set to piracy.whatever, piracy=1, public=0
    // - one with whatever set to public.whatever, piracy=0, public=1
    // The data might not, however, have any piracy and/or public data, in which case we omit the
    // relevant row.

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

    // Increase the matrix by 2520 row increments: 2520 is the lowest integer divisible by all
    // numbers from 1 to 10, so this will work properly for adding anywhere from 1 to 10 rows of
    // data per input row.  (Currently we only have 2 or 3, but 2520 is a reasonably sized number
    // for this, as well).
    rawdata.resize(2520, next_col);

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

    int rows = 0;
    for (const auto &row : csv) {
        if (rows >= rawdata.rows()) {
            rawdata.conservativeResize(rawdata.rows() + 2520, Eigen::NoChange);
        }
        // First row: pre. values
        for (unsigned fc = 0; fc < csv.fields().size(); fc++) {

            auto &d = csv.fields()[fc];
            if (d.substr(0, 7) != "piracy." and d.substr(0, 7) != "public.") {
                rawdata(rows, data_column[d]) = row[fc];
            }
        }
        if (piracy_data) rawdata(rows, data_column["piracy"]) = 0;
        if (public_data) rawdata(rows, data_column["public"]) = 0;
        for (auto &nan : pre_nan_) rawdata(rows, nan.second) = std::numeric_limits<double>::quiet_NaN();
        rows++;

        // Second row: piracy values
        if (piracy_data) {
            for (unsigned fc = 0; fc < csv.fields().size(); fc++) {
                auto &d = csv.fields()[fc];
                if (d.substr(0, 4) != "pre." and d.substr(0, 7) != "public.") {
                    rawdata(rows, data_column[d]) = row[fc];
                }
            }
            rawdata(rows, data_column["piracy"]) = 1;
            if (public_data) rawdata(rows, data_column["public"]) = 0;
            for (auto &nan : piracy_nan_) rawdata(rows, nan.second) = std::numeric_limits<double>::quiet_NaN();
            rows++;
        }
        // Third row: public values
        if (public_data) {
            for (unsigned fc = 0; fc < csv.fields().size(); fc++) {
                auto &d = csv.fields()[fc];
                if (d.substr(0, 4) != "pre." and d.substr(0, 7) != "piracy.") {
                    rawdata(rows, data_column[d]) = row[fc];
                }
            }
            if (piracy_data) rawdata(rows, data_column["piracy"]) = 0;
            rawdata(rows, data_column["public"]) = 1;
            for (auto &nan : public_nan_) rawdata(rows, nan.second) = std::numeric_limits<double>::quiet_NaN();
            rows++;
        }
    }

    // Resize the probably-oversized data matrix back to the number of rows we actually filled
    if (rawdata.rows() != rows) {
        rawdata.conservativeResize(rows, Eigen::NoChange);
    }

    for (auto &dc : data_column) {
        if (dc.first.substr(0, 4) == "pre." or dc.first.substr(0, 7) == "public." or dc.first.substr(0, 7) == "piracy.")
            continue;
        std::string name(dc.first);
        if (dc.first[0] == '.') name = dc.first.substr(1);
        data_.insert({name, SimpleVariable::create(name, rawdata.col(dc.second))});
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
                std::cerr << "Error: found NaN for " << avg_effects.equations()[i].depVar()->name() << ", row " <<
                    (r / (piracy_data and public_data ? 3 : piracy_data or public_data ? 2 : 1)) << "\n";
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
        for (auto &x : {"density", "cost_fixed", "cost_unit", "creation_time"}) {
            if (not data.has(x)) continue;
            if (not pre_nan_.count(x)) eq % data[x];
            if (piracy_data and not piracy_nan_.count(x)) eq % (data["piracy"] * data[x]);
            if (public_data and not public_nan_.count(x)) eq % (data["public"] * data[x]);
        }
        // These don't get included in "pre":
        for (auto &x : {"cost_piracy", "piracy_link_proportion"}) {
            if (not data.has(x)) continue;
            if (piracy_data and not piracy_nan_.count(x)) eq % (data["piracy"] * data[x]);
            if (public_data and not public_nan_.count(x)) eq % (data["public"] * data[x]);
        }
        // These don't get included in "pre" or "piracy":
        for (auto &x : {"public_sharing_tax"}) {
            if (not data.has(x)) continue;
            if (public_data and not public_nan_.count(x)) eq % (data["public"] * data[x]);
        }
        marg_effects.add(eq);
    }

    marg_effects.solve();
    std::cout << "\n\n\nMarginal effects:\n=================\n" << marg_effects;
}

