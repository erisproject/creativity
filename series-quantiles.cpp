#include "creativity/data/CSVParser.hpp"
#include "creativity/cmdargs/SeriesQuantiles.hpp"
#include "creativity/data/quantiles.hpp"
#include <boost/filesystem/operations.hpp>
#include <exception>
#include <string>
#include <iostream>

namespace creativity { namespace state { class State; } }

using namespace creativity;
using namespace creativity::data;
using namespace eris;

namespace fs = boost::filesystem;

int main(int argc, char *argv[]) {
    cmdargs::SeriesQuantiles args;
    try {
        args.parse(argc, argv);
    }
    catch (const std::exception &e) {
        std::cerr << "\n" << e.what() << "\n\n";
        exit(5);
    }

    std::istringstream iss(args.quantiles);
    std::set<double> quantiles;
    bool invalid = false;
    std::string quantile;
    while (std::getline(iss, quantile, ',')) {
        std::size_t pos;
        if (quantile == "max" or quantile == "maximum") quantiles.insert(1);
        else if (quantile == "min" or quantile == "minimum") quantiles.insert(0);
        else if (quantile == "median") quantiles.insert(0.5);
        else {
            try {
                double q = std::stod(quantile, &pos);
                if (pos != quantile.size() or q < 0 or q > 1) throw std::invalid_argument("invalid quantile value");
                quantiles.insert(q);
            } catch (const std::invalid_argument&) {
                std::cerr << "Error: requested quantile `" << quantile << "' is invalid\n";
                invalid = true;
            }
        }
    }

    if (invalid) {
        std::cerr << "Invalid quantile(s) provided; aborting.\n\n";
        exit(1);
    }
    else if (quantiles.empty()) {
        std::cerr << "No quantiles provided; aborting.\n\n";
        exit(1);
    }

    std::string output_header;
    {
        std::ostringstream headerss;
        headerss << "t";
        for (auto q : quantiles) headerss << "," << quantile_field(q);
        headerss << "\n";
        output_header = headerss.str();
    }

    if (args.output_prefix == args.output_unprefix) {
        std::cerr << "Invalid arguments: --prefix and --unprefix cannot be the same.\n\n";
        exit(4);
    }

    if (args.input.empty()) {
        std::cerr << "No series input files specified!\n\n";
        exit(2);
    }

    std::vector<std::string> output;
    output.reserve(args.input.size());
    for (const auto &input : args.input) {
        fs::path input_path(input);
        if (not fs::is_regular_file(input_path)) {
            std::cerr << "Error: input file `" << input << "' does not exist or is not a regular file\n";
            exit(3);
        }
        fs::path parent = input_path.parent_path();
        std::string filename = input_path.filename().string();
        bool changed_name = false;
        if (not args.output_unprefix.empty() and filename.substr(0, args.output_unprefix.size()) == args.output_unprefix) {
            filename = filename.substr(args.output_unprefix.size());
            changed_name = true;
        }
        if (not args.output_prefix.empty()) {
            filename = args.output_prefix + filename;
            changed_name = true;
        }
        if (not changed_name) {
            std::cerr << "Error: --prefix/--unprefix settings didn't imply a different output filename for `" << input << "': aborting.\n\n";
            exit(6);
        }
        fs::path output_path = parent/filename;
        std::string output_file = output_path.string();

        std::cout << input << " -> " << output_file << "..." << std::flush;

        CSVParser parser(input);
        std::ostringstream output_data;
        output_data << output_header;
        try {
            if (parser.fields().size() == 0) throw std::invalid_argument("no data fields found (not even t)");
            if (parser.fields()[0] != "t") throw std::invalid_argument("first field != t");
            if (parser.fields().size() < 2) throw std::invalid_argument("file contains no data (only a t column was found)");

            for (const auto &row : parser) {
                // Extract finite values (ignore any NaNs or infinities)
                std::vector<double> finite_values;
                std::copy_if(row.begin()+1, row.end(), std::back_inserter(finite_values), [](const double &v) { return std::isfinite(v); });
                if (finite_values.empty()) continue; // Completely skip time periods with zero finite values (often initial rows, possibly others)
                std::sort(finite_values.begin(), finite_values.end());

                output_data << row[0];
                for (auto q : quantiles) {
                    output_data << "," << double_str(data::quantile(finite_values.begin(), finite_values.end(), q), args.double_precision);
                }
                output_data << "\n";
            }
        }
        catch (const std::invalid_argument &e) {
            std::cerr << "\n\nError: input file `" << input << "' does not appear to be a valid creativity-series file: " << e.what() << ".  Aborting.\n\n";
            exit(8);
        }

        if (fs::exists(output_path)) {
            if (fs::equivalent(output_path, input_path)) {
                std::cerr << "\n\nError: output file and input file are the same for input file `" << input << "'; aborting.\n\n";
                exit(7);
            }
            else if (not args.overwrite) {
                std::cerr << "\n\nError: output file `" << output_file << "' already exists, and --overwrite was not specified.  Aborting.\n\n";
                exit(7);
            }
        }

        std::ofstream out(output_file, std::ios::out | std::ios::trunc);
        out << output_data.str();
        out.close();
        std::cout << " done.\n";
    }
}

