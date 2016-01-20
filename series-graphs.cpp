#include "creativity/data/CSVParser.hpp"
#include "creativity/cmdargs/SeriesGraphs.hpp"
#include "creativity/data/quantiles.hpp"
#include <eris/types.hpp>
#include <boost/filesystem/operations.hpp>
#include <cerrno>
#include <exception>
#include <string>
#include <iostream>

namespace creativity { namespace state { class State; } }

using namespace creativity;
using namespace creativity::data;
using namespace eris;

namespace fs = boost::filesystem;

int main(int argc, char *argv[]) {
    cmdargs::SeriesGraphs args;
    try {
        args.parse(argc, argv);
    }
    catch (const std::exception &e) {
        std::cerr << "\n" << e.what() << "\n\n";
        exit(5);
    }

    std::istringstream iss(args.levels);
    std::set<double> levels;
    bool invalid = false;
    std::string level;
    while (std::getline(iss, level, ',')) {
        std::size_t pos;
        if (level == "median") levels.insert(0);
        else {
            try {
                double q = std::stod(level, &pos);
                if (pos != level.size() or q < 0 or q > 1) throw std::invalid_argument("");
                levels.insert(q);
            } catch (const std::invalid_argument&) {
                std::cerr << "Error: requested level `" << level << "' is invalid\n";
                invalid = true;
            }
        }
    }

    if (invalid) {
        std::cerr << "Invalid confidence level(s) provided; aborting.\n\n";
        exit(1);
    }
    else if (levels.empty()) {
        std::cerr << "No confidence levels provided; aborting.\n\n";
        exit(1);
    }

    // Here's what we do now:
    // - Load each input file
    // - For series, calculate the needed quantiles
    // - Using the largest confidence level, calculate the minimum (from the lower quantile) and maximum (from the upper quantile).
    // - Find the lowest and highest "t" observations with data
    // - Repeat the above two steps across all input files
    // - This gives the x region (t range) and y region (data range) for the plots.
    // - For each file, plot the data:
    //   - In descending order of confidence levels, plot the ribbon associated with the confidence
    //     level; for 0 (the median) just plot a line instead.

    struct limits { double low, high; limits(double l, double h) : low{l}, high{h} {} };
    struct data {
        std::string filename;
        // Map t -> { level -> [lower, upper] }, with levels in descending order
        std::map<eris_time_t, std::map<double, limits, std::greater<double>>> quantiles;
        explicit data(std::string file) : filename{std::move(file)} {}
    };
    std::list<data> input;

    // DEBUG:
    auto trim = common_ends(args.input.begin(), args.input.end());
    std::cout << "Trimmed names:\n";
    for (const auto &file : args.input) {
        std::cout << "    " << file << ": (" << trim.first << "," << trim.second << ")" <<  file.substr(trim.first, file.length()-trim.second-trim.first) << "\n";
    }

    for (const auto &file : args.input) {
        fs::path input_path(file);
        if (not fs::is_regular_file(input_path)) {
            std::cerr << "Error: input file `" << file << "' does not exist or is not a regular file\n";
            exit(3);
        }

        CSVParser parser(file);
        input.emplace_back(file);
        auto &qdata = input.back().quantiles;

        try {
            if (parser.fields()[0] != "t") throw std::invalid_argument("first field != t");

            // First figure out whether this is a series or quantiles file:
            bool found_series = false, found_quantiles = false;
            for (size_t i = 1; i < parser.fields().size(); i++) {
                if (not found_quantiles and std::regex_match(parser.fields()[i], ordinal_regex))
                    found_series = true;
                else if (not found_series and std::regex_match(parser.fields()[i], quantile_field_regex))
                    found_quantiles = true;
                else {
                    found_series = found_quantiles = false;
                    break;
                }
            }
            if (not found_series and not found_quantiles)
                throw std::invalid_argument("file headers do not match either a series or quantiles file");

            if (found_quantiles) {
                // Map quantile to field index:
                std::unordered_map<double, size_t> l_pos_low, l_pos_high;

                // Make sure we have all the quantiles we need
                for (const auto &l : levels) {
                    double qlow = 0.5 - l/2, qhigh = 0.5 + l/2;
                    for (double q : {qlow, qhigh}) {
                        auto qheader = quantile_field(q);
                        if (not parser.hasField(qheader))
                            throw std::invalid_argument("file is missing required " + qheader + " quantile field");
                    }
                    l_pos_low.emplace (l, parser.fieldPosition(quantile_field(qlow)));
                    l_pos_high.emplace(l, parser.fieldPosition(quantile_field(qhigh)));
                }

                for (const auto &row : parser) {
                    for (const auto &l : levels) {
                        qdata[row[0]].emplace(l, limits(row[l_pos_low[l]], row[l_pos_high[l]]));
                    }
                }
            }
            else { // found_series
                parser.allow_missing_values = parser.header().size()-1;

                for (const auto &row : parser) {
                    if (row.size() == 1) continue; // If we have just the time index but no data at all for that index skip it.

                    // Make sure the data is sorted:
                    for (int i = 2; i < row.size(); i++) {
                        if (row[i] < row[i-1]) throw std::invalid_argument("file contains unsorted data (line "
                                + std::to_string(parser.lineNumber()) + ", columns " + std::to_string(i-1) + "--" + std::to_string(i) + ")");

                    }

                    for (const auto &l : levels) {
                        qdata[row[0]].emplace(l, limits(
                                    creativity::data::quantile(row.tail(row.size()-1), 0.5 - l/2),
                                    creativity::data::quantile(row.tail(row.size()-1), 0.5 + l/2)));
                    }
                }
            }
        }
        catch (const std::invalid_argument &e) {
            std::cerr << "\n\nError: input file `" << file << "' does not appear to be a valid series or quantiles file: " << e.what() << ".  Aborting.\n\n";
            exit(8);
        }
    }

    double ymin = std::numeric_limits<double>::infinity(), ymax = -std::numeric_limits<double>::infinity();

    eris_time_t tmin = std::numeric_limits<eris_time_t>::max(), tmax = std::numeric_limits<eris_time_t>::min();
    for (auto &f : input) for (auto &period : f.quantiles) {
        if (period.first < tmin) tmin = period.first;
        if (period.first > tmax) tmax = period.first;

        auto &lims = period.second.begin()->second;
        if (lims.low < ymin) ymin = lims.low;
        if (lims.high > ymax) ymax = lims.high;
    }

    std::cout << "Found t range: [" << tmin << ", " << tmax << "], y range: [" << ymin << ", " << ymax << "]\n";
}

