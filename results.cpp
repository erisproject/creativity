#include "creativity/data/CSVParser.hpp"
#include "creativity/data/SUR.hpp"
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
std::unordered_map<std::string, SimpleVariable> data;

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
    bool piracy_rows = false, public_rows = false;

    int next_col = 0;
    for (auto &f : csv.fields()) {
        if (f[0] == '.') throw std::runtime_error("CSV file has invalid field `" + f + "': fields may not start with .");
        std::string data_field;
        if (f.substr(0, 4) == "pre.") data_field = f.substr(3);
        else if (f.substr(0, 7) == "piracy.") { data_field = f.substr(6); piracy_rows = true; }
        else if (f.substr(0, 7) == "public.") { data_field = f.substr(6); public_rows = true; }

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

    if (piracy_rows) data_column.insert({"piracy", next_col++}); // piracy dummy
    if (public_rows) data_column.insert({"public", next_col++}); // public sharing dummy

    // Increase the matrix by 2520 row increments: 2520 is the lowest integer divisible by all
    // numbers from 1 to 10, so this will work properly for adding anywhere from 1 to 10 rows of
    // data per input row.  (Currently we only have 2 or 3, but 2520 is a reasonably sized number
    // for this, as well).
    rawdata.resize(2520, next_col);

    // Some fields don't have all three (for example, books_pirated has piracy_ and public_ values
    // but not a pre_ value; books_public_copies is only under public_), so we need to fill in some
    // NaNs for the ones where they don't exist: track the columns needing NaNs here:
    std::list<int> pre_nans, piracy_nans, public_nans;
    for (auto &pair : data_column) {
        if (pair.first[0] == '.') {
            if (data_column.count("pre" + pair.first) == 0)
                pre_nans.push_back(pair.second);
            if (piracy_rows and data_column.count("piracy" + pair.first) == 0)
                piracy_nans.push_back(pair.second);
            if (public_rows and data_column.count("public" + pair.first) == 0)
                public_nans.push_back(pair.second);
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
        if (piracy_rows) rawdata(rows, data_column["piracy"]) = 0;
        if (public_rows) rawdata(rows, data_column["public"]) = 0;
        for (auto nan_i : pre_nans) rawdata(rows, nan_i) = std::numeric_limits<double>::quiet_NaN();
        std::cerr << "after row, rawdata row: " << rawdata.row(rows) << "\n";
        rows++;

        // Second row: piracy values
        if (piracy_rows) {
            for (unsigned fc = 0; fc < csv.fields().size(); fc++) {
                auto &d = csv.fields()[fc];
                if (d.substr(0, 4) != "pre." and d.substr(0, 7) != "public.") {
                    rawdata(rows, data_column[d]) = row[fc];
                }
            }
            rawdata(rows, data_column["piracy"]) = 1;
            if (public_rows) rawdata(rows, data_column["public"]) = 0;
            for (auto nan_i : piracy_nans) rawdata(rows, nan_i) = std::numeric_limits<double>::quiet_NaN();
            rows++;
        }
        // Third row: public values
        if (public_rows) {
            for (unsigned fc = 0; fc < csv.fields().size(); fc++) {
                auto &d = csv.fields()[fc];
                if (d.substr(0, 4) != "pre." and d.substr(0, 7) != "piracy.") {
                    rawdata(rows, data_column[d]) = row[fc];
                }
            }
            if (piracy_rows) rawdata(rows, data_column["piracy"]) = 0;
            rawdata(rows, data_column["public"]) = 1;
            for (auto nan_i : public_nans) rawdata(rows, nan_i) = std::numeric_limits<double>::quiet_NaN();
            rows++;
        }
    }

    // Resize the probably-oversized data matrix back to the number of rows we actually filled
    if (rawdata.rows() != rows) {
        rawdata.conservativeResize(rows, Eigen::NoChange);
    }

    std::cout << "rawdata:\n" << rawdata << "\n";



    for (auto &dc : data_column) {
        if (dc.first.substr(0, 4) == "pre." or dc.first.substr(0, 7) == "public." or dc.first.substr(0, 7) == "piracy.")
            continue;
        std::string name(dc.first);
        if (dc.first[0] == '.') name = dc.first.substr(1);
//        SimpleVariable var(name, rawdata, dc.second);
        data.insert({name, SimpleVariable(name, rawdata, dc.second)});
    }

    for (auto &dc : data) {
        std::cout << dc.first << ":\n" << dc.second.values() << "\n\n";
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

    // TODO: it would be nice to be able to specify these models via cli args instead of hard-coding
    // the model here
//    SUR avg_piracy(
//            Equation(data["utility"]) + 1 + 
}

