#include "creativity/Creativity.hpp"
#include "creativity/data/simdata.hpp"
#include "creativity/data/util.hpp"
#include "creativity/state/Storage.hpp"
#include "creativity/cmdargs/Series.hpp"
#include <eris/types.hpp>
#include <boost/filesystem/operations.hpp>
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
#include <vector>
#include <set>
#include <fstream>
#include <sstream>

namespace creativity { namespace state { class State; } }

using namespace creativity;
using namespace creativity::state;
using namespace eris;

// The series we want to calculate, as given to --series
std::unordered_set<std::string> series_wanted;
// Results: "books_written" -> array of periods -> array of values
std::unordered_map<std::string, std::vector<std::multiset<double>>> values;
unsigned int values_count, error_count;
std::mutex values_mutex; // Guards values, output_count, error_count
decltype(cmdargs::Series::input.cbegin()) input_it, input_it_end;
std::mutex input_it_mutex;

// These have to agree across files:
eris_time_t periods = 0, piracy_begins = 0, public_begins = 0;
bool allow_unused_periods = false; // Will only be true if --periods is explicitly given
bool need_parameters = true;


// Parses a file; takes the program arguments and a vector of datum objects used to calculate the
// values of interest, which are stored in `values`.
void thr_parse_file(
        const cmdargs::Series &args,
        const std::vector<data::datum> &data) {

    input_it_mutex.lock();

    while ((error_count == 0 or args.ignore_errors) and input_it != input_it_end) {
        std::string  source(*input_it++);
        // Hold the lock if this is the first file and we weren't given all of the needed simulation
        // period values: we have to set periods, piracy_begins, public_begins (the lock also
        // protects assigning to these variables).
        if (not need_parameters) input_it_mutex.unlock();

        values_mutex.lock();
        std::cout << "Processing " << source << "\n";
        values_mutex.unlock();

#define FAIL(WHY) \
        values_mutex.lock(); \
        std::cerr << "Error reading `" << source << "': " << WHY << "\n"; \
        error_count++; \
        values_mutex.unlock(); \
        continue;

        std::ostringstream output;
        output.precision(args.double_precision);
        auto creativity = Creativity::create();
        // Filename input
        try {
            creativity->fileRead(source);
        }
        catch (std::ios_base::failure&) {
            FAIL("I/O error: " << std::strerror(errno));
        }
        catch (std::exception &e) {
            FAIL("An exception occured: " << e.what());
        }

        auto locked_storage = creativity->storage();
        auto &storage = *locked_storage.first;

        // If parameters were given, use them; otherwise infer them from the first file read
        if (need_parameters) {
            if (args.periods != 0) {
                periods = args.periods;
                allow_unused_periods = args.allow_unused_periods;
            }
            else {
                periods = storage.size()-1;
                std::cout << "Inferred --periods " << periods << "\n";
            }
            if (args.piracy_begins != 0) piracy_begins = args.piracy_begins;
            else {
                piracy_begins = creativity->parameters.piracy_begins;
                std::cout << "Inferred --piracy-begins " << piracy_begins << "\n";
            }
            if (args.public_sharing_begins != 0) public_begins = args.public_sharing_begins;
            else {
                public_begins = creativity->parameters.public_sharing_begins;
                std::cout << "Inferred --public-sharing-begins " << public_begins << "\n";
            }
            need_parameters = false;
            input_it_mutex.unlock();
        }

        // If we need more than is available, we can't use this file.
        if (allow_unused_periods ? periods+1 > storage.size() : periods+1 != storage.size()) {
            FAIL(periods << " periods required, but file has " << storage.size()-1);
        }
        // Piracy begins doesn't have to match if the requested number of periods is less than the
        // piracy begins value, but otherwise does:
        if (piracy_begins <= periods and piracy_begins != creativity->parameters.piracy_begins) {
            FAIL("simulation piracy begins at t=" << creativity->parameters.piracy_begins << " but t=" << piracy_begins << " is required");
        }
        // Likewise for public sharing
        if (public_begins <= periods and public_begins != creativity->parameters.public_sharing_begins) {
            FAIL("simulation public sharing begins at t=" << creativity->parameters.public_sharing_begins << " but t=" << public_begins << " is required");
        }

        // Calculate all the local values, then copy them into the results variable in one shot
        std::unordered_map<std::string, std::vector<double>> local_values;
        for (eris_time_t t = 0; t <= periods; t++) {
            // Cache the period values so it stays in memory and doesn't have to be re-read each
            // time
            std::shared_ptr<const State> cached = storage[t];

            for (const auto &datum : data) {
                if (series_wanted.count(datum.name) > 0) {
                    if (t == 0) local_values[datum.name].resize(periods+1, std::numeric_limits<double>::quiet_NaN());
                    local_values[datum.name][t] = datum.calculate(storage, t, t);
                }
            }
        }

        // Copy values we read into the overall values container
        values_mutex.lock();
        for (const auto &v : local_values) {
            auto &vstore = values[v.first];
            if (vstore.size() < v.second.size()) vstore.resize(v.second.size());
            for (unsigned t = 0; t < v.second.size(); t++) {
                // Only copy finite values (this mainly excluded NaNs as nothing should be capable
                // of generating infinity).
                double value = v.second[t];
                if (std::isfinite(value)) vstore[t].insert(value);
            }
        }
        values_count++;
        values_mutex.unlock();

        input_it_mutex.lock();
    }
    input_it_mutex.unlock();
}

int main(int argc, char *argv[]) {
    cmdargs::Series args;
    try {
        args.parse(argc, argv);
    }
    catch (const std::exception &e) {
        std::cerr << "\n" << e.what() << "\n\n";
        exit(5);
    }

    // Simulation fields:
    auto data = data::data_fields();

    if (args.help_series) {
        std::cout << args.version() << "\n\n";
        std::cout << "Supported series parameters:\n\n";

        for (auto &d : data) {
            std::cout << "        " << d.name << "\n";
        }

        std::cout << "\n\nNote that some values (e.g. book_quality) can be NaN in certain circumstances\n"
            "(for example, no books written during a period): any such NaN values are excluded\n"
            "from the generated data.  In effect this means that the length of the list of\n"
            "values for each period can change.\n";

        exit(0);
    }

    // Verify the series requested
    {
        std::unordered_set<std::string> series_available;
        for (auto &d : data) {
            series_available.emplace(d.name);
        }

        std::istringstream iss(args.series);
        std::string series;
        bool invalid = false;
        while (std::getline(iss, series, ',')) {
            if (series_available.count(series)) {
                series_wanted.insert(series);
            }
            else {
                std::cerr << "Error: requested variable series `" << series << "' is not a valid simulation variable\n";
                invalid = true;
            }
        }

        if (not invalid and series_wanted.empty()) {
            std::cerr << "Error: no variable series specified\n";
            invalid = true;
        }

        if (invalid) exit(7);
    }


    try {
        boost::filesystem::create_directories(args.output_dir);
    }
    catch (const std::exception &e) {
        std::cerr << "\nUnable to create output directory: " << e.what() << "\n\n";
        exit(1);
    }

    input_it = args.input.cbegin();
    input_it_end = args.input.cend();

    if (input_it == input_it_end) {
        std::cerr << "\nError: no simulation input files specified\n\n";
        exit(1);
    }

    if (args.threads == 0) {
        thr_parse_file(args, data);
    }
    else {
        Eigen::initParallel();
        std::vector<std::thread> threads;
        // Start up the requested number of threads:
        for (unsigned t = 0; t < args.threads; t++) {
            threads.emplace_back(thr_parse_file, args, data);
        }
        // Wait for all threads to finish:
        for (auto &th : threads) {
            if (th.joinable()) th.join();
        }
    }

    if (error_count > 0) {
        if (args.ignore_errors) {
            std::cerr << "Warning: some files failed to be read or were unsuitable.\n";
        }
        else {
            std::cerr << "Error: encountered unusable file; aborting.\n";
            exit(2);
        }
    }

    if (values_count == 0) {
        std::cerr << "Error: no usable data sources!\n";
        exit(1);
    }

    std::cout << "Successfully processed " << values_count << "/" << (values_count+error_count) << " simulation files.\n\n";

    for (auto &v : values) {
        // Figure out how many "nth" headers we need for this series:
        unsigned n = 0;
        for (auto &time : v.second) {
            if (time.size() > n) n = time.size();
        }

        // Write an output file
        std::string output_file = args.output_dir + "/series-" + v.first + ".csv";
        std::ofstream out(output_file, std::ios::out | std::ios::trunc);
        out << "t";
        for (unsigned i = 1; i <= n; i++) {
            out << "," << i << (
                (i % 10 == 1 and not i % 100 == 11) ? "st" :
                (i % 10 == 2 and not i % 100 == 12) ? "nd" :
                (i % 10 == 3 and not i % 100 == 13) ? "rd" :
                "th");
        }
        out << "\n";
        for (eris_time_t t = 0; t < v.second.size(); t++) {
            out << t;
            for (const auto &val : v.second[t]) {
                out << "," << data::double_str(val, args.double_precision);
            }
            out << "\n";
        }

        out.close();
        std::cout << "Wrote '" << v.first << "' series data to " << output_file << "\n";
    }
}

