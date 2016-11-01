#include "creativity/Creativity.hpp"
#include "creativity/data/simdata.hpp"
#include "creativity/data/util.hpp"
#include "creativity/state/Storage.hpp"
#include "creativity/state/FileStorage.hpp"
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

using Locker = std::lock_guard<std::mutex>;

// The series we want to calculate, as given to --series
std::unordered_set<std::string> series_wanted;

// Values container.  values["field_name"][t][fileindex] = value
std::unordered_map<std::string, std::vector<std::vector<double>>> values;
// File list, in same order as [fileindex] above:
std::vector<std::string> files;
size_t error_count = 0;
std::mutex values_mutex; // Guards values, files, error_count

std::mutex input_mutex; // Guards all of the below:
size_t input_index;
std::queue<std::pair<std::string, std::unique_ptr<std::stringstream>>> preload_queue;
bool preload_done = false;
std::condition_variable preload_cv; // Used to wait for preload data to appear
std::condition_variable preload_next_cv; // Used to tell the preload thread that something has been removed from the preload queue
// These have to agree across files:
eris_time_t periods = 0, piracy_begins = 0, policy_begins = 0;
bool allow_unused_periods = false; // Will only be true if --periods is explicitly given
bool need_parameters = true;

void thr_preload(const cmdargs::Series &args) {
    std::unique_lock<std::mutex> input_lock(input_mutex);
    const unsigned int max_queue = args.preload;
    std::cerr << "preload = " << max_queue << "\n";
    if (max_queue < 1) throw std::logic_error("Internal error: thr_preload called in non-preload mode");
    while ((error_count == 0 or args.ignore_errors) and input_index < args.input.size()) {
        preload_next_cv.wait(input_lock, [&max_queue]{ return preload_queue.size() < max_queue; });
        while (preload_queue.size() < max_queue && input_index < args.input.size()) {
            std::string source(args.input[input_index++]);
            input_lock.unlock();
            std::unique_ptr<std::stringstream> s(new std::stringstream(
                        std::ios_base::in | std::ios_base::out | std::ios_base::binary));
            s->exceptions(s->failbit | s->badbit);
            bool success = false;
            std::ifstream f;
            try {
                f.exceptions(f.badbit | f.failbit);
                f.open(source, std::ios_base::in | std::ios_base::binary);
                *s << f.rdbuf();
                success = true;
            }
            catch (std::ios_base::failure &e) {
                std::lock_guard<std::mutex> g(values_mutex);
                std::cerr << "Unable to preload `" << source << "': " << std::strerror(errno) << "\n";
                error_count++;
            }
            catch (std::exception &e) {
                std::lock_guard<std::mutex> g(values_mutex);
                std::cerr << "Unable to preload `" << source << "': " << e.what() << "\n";
                error_count++;
            }

            input_lock.lock();
            if (success) {
                std::cout << "preloaded " << source << "\n";
                preload_queue.emplace(std::move(source), std::move(s));
                preload_cv.notify_all();
            }
        }
    }

    std::cout << "preload done\n";

    // Input file preloading is done
    preload_done = true;
    preload_cv.notify_all();
    input_lock.unlock();
}

// Parses a file; takes the program arguments and a vector of datum objects used to calculate the
// values of interest, which are stored in `values`.
void thr_parse_file(
        const cmdargs::Series &args,
        const std::vector<data::datum> &data) {

    std::unique_lock<std::mutex> input_lock(input_mutex);
    const bool preload_mode = args.preload > 0;
    size_t input_i;
    while ((error_count == 0 or args.ignore_errors) and (preload_mode || (input_i = input_index++) < args.input.size())) {
        std::string source;
        std::unique_ptr<std::stringstream> ss;
        if (preload_mode) {
            preload_cv.wait(input_lock, []{ return preload_done || !preload_queue.empty(); });
            if (!preload_queue.empty()) {
                {
                    auto &next = preload_queue.front();
                    source = std::move(next.first);
                    ss = std::move(next.second);
                }
                preload_queue.pop();
                preload_next_cv.notify_all();
                if (!ss) {
                    std::lock_guard<std::mutex> g(values_mutex);
                    std::cerr << "Unable to read preloaded file `" << source << "': stringstream pointer is null\n";
                    continue;
                }
            }
            else {
                // Preloading is finished *and* the queue is empty, so we're done.
                break;
            }
        }
        else {
            source = args.input[input_i];
        }

        // Hold the lock if this is the first file and we weren't given all of the needed simulation
        // period values: we have to set periods, piracy_begins, policy_begins (the lock also
        // protects assigning to these variables) before unlocking; this essentially means that only
        // one thread runs until the first thread has determined that initial info.
        if (not need_parameters) input_lock.unlock();

        { Locker l(values_mutex); std::cout << "Processing " << source << std::endl; }

#define FAIL(WHY) { Locker l(values_mutex); std::cerr << "Error reading `" << source << "': " << WHY << std::endl; error_count++; continue; }

        Creativity creativity;
        // Filename input
        try {
            if (preload_mode) // Preloaded input:
                creativity.read<FileStorage>(std::move(ss), FileStorage::Mode::READONLY);
            else
                creativity.read<FileStorage>(source, FileStorage::Mode::READONLY, args.memory_xz, args.tmpdir);
        }
        catch (std::ios_base::failure&) FAIL("I/O error: " << std::strerror(errno))
        catch (std::exception &e) FAIL("An error occured: " << e.what())

        auto locked_storage = creativity.storage();
        auto &storage = *locked_storage.first;

        // If parameters were given, use them; otherwise infer them from the first file read
        if (need_parameters) {
            if (args.periods != 0) {
                periods = args.periods;
                allow_unused_periods = args.allow_unused_periods;
            }
            else {
                periods = storage.size()-1;
                std::cout << "Inferred --periods " << periods << std::endl;
            }
            if (args.piracy_begins != 0) piracy_begins = args.piracy_begins;
            else {
                piracy_begins = creativity.parameters.piracy_begins;
                std::cout << "Inferred --piracy-begins " << piracy_begins << std::endl;
            }
            if (args.policy_begins != 0) policy_begins = args.policy_begins;
            else {
                policy_begins = creativity.parameters.policy_begins;
                std::cout << "Inferred --policy-begins " << policy_begins << std::endl;
            }
            need_parameters = false;
            input_lock.unlock();
        }

        // If we need more than is available, we can't use this file.
        if (allow_unused_periods ? periods+1 > storage.size() : periods+1 != storage.size()) {
            FAIL(periods << " periods required, but file has " << storage.size()-1);
        }
        // Piracy begins doesn't have to match if the requested number of periods is less than the
        // piracy begins value, but otherwise does:
        if (piracy_begins <= periods and piracy_begins != creativity.parameters.piracy_begins) {
            FAIL("simulation piracy begins at t=" << creativity.parameters.piracy_begins << " but t=" << piracy_begins << " is required");
        }
        // Likewise for policy
        if (policy_begins <= periods and policy_begins != creativity.parameters.policy_begins) {
            FAIL("simulation policy begins at t=" << creativity.parameters.policy_begins << " but t=" << policy_begins << " is required");
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
        {
            Locker locker(values_mutex);
            auto fileindex = files.size();
            files.push_back(source);

            // values["field_name"][t][fileindex] = value
            for (const auto &v : local_values) {
                auto &vstore = values[v.first];
                if (vstore.size() < v.second.size()) vstore.resize(v.second.size());
                for (unsigned t = 0; t < v.second.size(); t++) {
                    vstore[t].reserve(args.input.size());
                    if (vstore[t].size() <= fileindex) vstore[t].resize(fileindex+1, std::numeric_limits<double>::quiet_NaN());
                    vstore[t][fileindex] = v.second[t];
                }
            }
        }

        input_lock.lock();
    }
    input_lock.unlock();
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

    if (args.input.empty()) {
        std::cerr << "\nError: no simulation input files specified\n\n";
        exit(50);
    }

    // Check for files specified more than once, and abort if found
    {
        std::unordered_set<std::string> seen;
        for (const auto &source : args.input) {
            auto ins = seen.insert(source);
            if (not ins.second) {
                std::cerr << "\nError: simulation input files contains duplicate entry `" << source << "'; aborting\n\n";
                exit(51);
            }
        }
    }

    files.reserve(args.input.size());

    std::vector<std::thread> threads;
    if (args.preload > 0) {
        // Start the preload thread, which loads files into memory:
        threads.emplace_back(thr_preload, args);
    }

    if (args.threads == 0) {
        thr_parse_file(args, data);
    }
    else {
        Eigen::initParallel();
        // Start up the requested number of threads:
        for (unsigned t = 0; t < args.threads; t++) {
            threads.emplace_back(thr_parse_file, args, data);
        }
    }

    // Wait for all threads to finish:
    for (auto &th : threads) {
        if (th.joinable()) th.join();
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

    if (files.empty()) {
        std::cerr << "Error: no usable data sources!\n";
        exit(1);
    }

    std::cout << "Successfully processed " << files.size() << "/" << (files.size()+error_count) << " simulation files.\n" << std::endl;

    std::string header;
    {
        // Maps actual source names into CSV-compatible source names by removing potentially problematic
        // characters, and appending a unique number to resolve any duplicates.
        std::unordered_set<std::string> source_used;
        std::ostringstream headeross;
        headeross << "t";
        for (const auto &source : files) {
            std::string fixed = data::csv_fix(source);
            std::string try_s = fixed;
            int append = 2;
            while (not source_used.insert(try_s).second) {
                try_s = fixed + "-" + std::to_string(append++);
            }
            headeross << "," << try_s;
        }
        headeross << "\n";
        header = headeross.str();
    }

    for (auto &v : values) {
        // Write an output file
        std::string output_file = args.output_dir + "/series-" + v.first + ".csv";
        std::ofstream out(output_file, std::ios::out | std::ios::trunc);
        out.exceptions(out.failbit | out.badbit);
        out << header;

        for (eris_time_t t = 0; t < v.second.size(); t++) {
            out << t;
            for (const auto &val : v.second[t]) {
                out << "," << data::double_str(val, args.double_precision);
            }
            out << "\n";
        }

        out.close();
        std::cout << "Wrote '" << v.first << "' series data to " << output_file << std::endl;
    }
}

