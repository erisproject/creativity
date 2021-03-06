#include "creativity/Creativity.hpp"
#include "creativity/data/simdata.hpp"
#include "creativity/state/Storage.hpp"
#include "creativity/state/FileStorage.hpp"
#include "creativity/cmdargs/Data.hpp"
#include <eris/types.hpp>
#include <cerrno>
#include <exception>
#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <list>
#include <queue>
#include <mutex>
#ifdef __linux__
extern "C" {
#include <sys/prctl.h>
}
#endif

namespace creativity { namespace state { class State; } }

using namespace creativity;
using namespace creativity::state;
using namespace eris;

// When outputting CSV values, we want exact precision, but don't necessarily need the full
// max_digits10 digits to get there.  For example, 0.1 with max_digits10 becomes 0.10000000000000001
// but 0.1 also converts to the numerically identical value.  0.10000000000000002, on the other
// hand, needs every decimal digit of precision, because it isn't equal to the double value 0.1
//
// This function produces tries first at the requested precision, then the requested precision less
// one, then less two; if the subsequent values are numerically identical to the given double value,
// the shortest is returned.
std::string double_str(double d, unsigned precision) {
    if (std::isfinite(d)) {
        double round_trip;
        for (unsigned prec : {precision-2, precision-1}) {
            std::stringstream ss;
            ss.precision(prec);
            ss << d;
            ss >> round_trip;
            if (round_trip == d) { ss << d; return ss.str(); }
        }
    }
    std::stringstream ss;
    ss.precision(precision);
    ss << d;
    return ss.str();
}

unsigned int processing_counter = 0, output_count = 0, input_count = 0;
std::mutex output_mutex; // Guards std::cout, std::cerr, and the above counters
decltype(cmdargs::Data::input.cbegin()) input_it, input_it_end;
std::queue<std::pair<std::string, std::unique_ptr<std::stringstream>>> preload_queue;
bool preload_done = false;
std::mutex input_mutex; // Guards both the above iterator and preload_queue
std::condition_variable preload_cv; // Used to wait for preload data to appear
std::condition_variable preload_next_cv; // Used to tell the preload thread that something has been removed from the preload queue

#ifdef __linux__
// Machinery for helping the main process managing the threads, so that we can let the parent update
// its status via prctl.
std::mutex thread_mutex;
unsigned threads_running = 0;
std::condition_variable thread_done_cv; // Signals the parent that a thread is done
#define THREAD_START { std::lock_guard<std::mutex> g(thread_mutex); threads_running++; }
#define THREAD_DONE { std::lock_guard<std::mutex> g(thread_mutex); threads_running--; thread_done_cv.notify_all(); }
#else
#define THREAD_START
#define THREAD_DONE
#endif

void thr_preload(const cmdargs::Data &args) {
    std::unique_lock<std::mutex> input_lock(input_mutex);
    const unsigned int max_queue = args.preload;
    if (max_queue < 1) throw std::logic_error("Internal error: thr_preload called in non-preload mode");
    while (input_it != input_it_end) {
        preload_next_cv.wait(input_lock, [&max_queue]{ return preload_queue.size() < max_queue; });
        while (preload_queue.size() < max_queue && input_it != input_it_end) {
            std::string source(*input_it++);
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
                std::lock_guard<std::mutex> g(output_mutex);
                std::cerr << "Unable to preload `" << source << "': " << std::strerror(errno) << "\n";
            }
            catch (std::exception &e) {
                std::lock_guard<std::mutex> g(output_mutex);
                std::cerr << "Unable to preload `" << source << "': " << e.what() << "\n";
            }

            input_lock.lock();
            if (success) {
                preload_queue.emplace(std::move(source), std::move(s));
                preload_cv.notify_all();
            }
        }
    }

    // Input file preloading is done
    preload_done = true;
    preload_cv.notify_all();
    input_lock.unlock();

    THREAD_DONE
}

void thr_parse_file(
        const cmdargs::Data &args,
        const std::vector<data::initial_datum> &initial_data,
        const std::vector<data::datum> &data,
        unsigned readable_name_width) {
    const bool preload_mode = args.preload > 0;
    std::unique_lock<std::mutex> input_lock(input_mutex);
    while (preload_mode || input_it != input_it_end) {
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
                    std::lock_guard<std::mutex> g(output_mutex);
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
            source = *input_it++;
        }
        input_lock.unlock();

        output_mutex.lock();
        processing_counter++;
        std::cerr << "Processing [" << processing_counter << "/" << input_count << "]: " << source << "\n";

        output_mutex.unlock();

        std::ostringstream output;
        output.precision(args.double_precision);
        Creativity creativity;
        try {
            if (preload_mode) // Preloaded input:
                creativity.read<FileStorage>(std::move(ss), FileStorage::Mode::READONLY);
            else // File input:
                creativity.read<FileStorage>(source, FileStorage::Mode::READONLY, args.memory_xz, args.tmpdir);
        }
        catch (std::ios_base::failure&) {
            std::lock_guard<std::mutex> g(output_mutex);
            std::cerr << "Unable to read/parse `" << source << "': " << std::strerror(errno) << "\n";
            continue;
        }
        catch (std::exception &e) {
            std::lock_guard<std::mutex> g(output_mutex);
            std::cerr << "Unable to read/parse `" << source << "': " << e.what() << "\n";
            continue;
        }

        auto locked_storage = creativity.storage();
        auto &storage = *locked_storage.first;

        // Check that the data source has enough data:

        const eris_time_t &piracy_begins = creativity.parameters.piracy_begins,
              &policy_begins = creativity.parameters.policy_begins;
        eris_time_t post_piracy = 0, post_policy = 0, post_pre = 0;
        if (not args.skip.piracy)
            post_piracy = (args.skip.policy ? storage.size() : policy_begins);
        if (not args.skip.policy)
            post_policy = storage.size();
        post_pre = std::min(storage.size(),
                args.skip.piracy
                ? args.skip.policy
                    ? storage.size()
                    : policy_begins
                : piracy_begins);

#define SKIP_IF(CONDITION, REASON) if (CONDITION) { \
    std::lock_guard<std::mutex> g(output_mutex); \
    std::cerr << "Skipping `" << source << "': " << REASON << "\n"; \
    continue; \
}

        SKIP_IF(args.verify.periods > 0 and (size_t) args.verify.periods != storage.size()-1,
                "simulation periods " << storage.size()-1 << " != " << args.verify.periods);
        SKIP_IF(args.verify.piracy_begins > 0 and (size_t) args.verify.piracy_begins != piracy_begins,
                "simulation piracy_begins " << piracy_begins << " != " << args.verify.piracy_begins);
        SKIP_IF(args.verify.policy_begins > 0 and (size_t) args.verify.policy_begins != policy_begins,
                "simulation policy_begins " << policy_begins << " != " << args.verify.policy_begins);

        SKIP_IF(post_pre-1 <= args.data_periods,
                "simulation \"pre\" periods " << post_pre-2 << " < " << args.data_periods);

        // If we're including both short-run and long-run effects, we need twice as many periods to
        // avoid overlap.
        unsigned periods_needed = args.data_periods;
        if (not args.skip.short_run) periods_needed *= 2;

        if (not args.skip.piracy) {
            SKIP_IF(post_piracy < piracy_begins + periods_needed,
                    "simulation doesn't have enough piracy periods (t=" << piracy_begins << " through t=" << post_piracy-1 << ")");
        }
        if (not args.skip.policy) {
            SKIP_IF(post_policy < policy_begins + periods_needed,
                    "simulation doesn't have enough policy periods (t=" << policy_begins << " through t=" << post_policy-1 << ")");
        }
#undef SKIP_IF

        if (args.human_readable) output << "\n\n" << source << "\n==========" << std::endl;
        else output << data::csv_fix(source);
        for (auto &d : initial_data) {
            if (args.human_readable) output << std::setw(readable_name_width+1) << d.name + ":" << " ";
            else output << ",";

            if (d.calc_double) output << double_str(d.calc_double(creativity.parameters), args.double_precision);
            else output << d.calc_int(creativity.parameters);

            if (args.human_readable) output << "\n";
        }

        // Keep a reference to all of the periods, so that the following don't end up continually
        // reloading them from disk/database.
        std::list<std::shared_ptr<const State>> state_cache;

        for (eris_time_t t = post_pre    - args.data_periods; t < post_pre;    t++) state_cache.push_back(storage[t]);

        // pre_*:
        for (auto &d : data) {
            if (d.applies_to.pre) {
                if (args.human_readable) output << std::setw(readable_name_width+1) << "pre." + d.name + ":" << " ";
                else output << ",";

                output << double_str(d.calculate(storage, post_pre - args.data_periods, post_pre - 1), args.double_precision);

                if (args.human_readable) output << "\n";
            }
        }


        // piracy.SR.*:
        // piracy.*:
        if (not args.skip.piracy) {
            if (not args.skip.short_run) {
                state_cache.clear();
                for (eris_time_t t = 0; t < args.data_periods; t++) state_cache.push_back(storage[piracy_begins + t]);

                for (auto &d : data) {
                    if (d.applies_to.piracy) {
                        if (args.human_readable) output << std::setw(readable_name_width+1) << "piracy.SR." + d.name + ":" << " ";
                        else output << ",";

                        output << double_str(d.calculate(storage, piracy_begins, piracy_begins + args.data_periods - 1), args.double_precision);

                        if (args.human_readable) output << "\n";
                    }
                }
            }

            state_cache.clear();
            for (eris_time_t t = post_piracy - args.data_periods; t < post_piracy; t++) state_cache.push_back(storage[t]);

            for (auto &d : data) {
                if (d.applies_to.piracy) {
                    if (args.human_readable) output << std::setw(readable_name_width+1) << "piracy." + d.name + ":" << " ";
                    else output << ",";

                    output << double_str(d.calculate(storage, post_piracy - args.data_periods, post_piracy - 1), args.double_precision);

                    if (args.human_readable) output << "\n";
                }
            }
        }

        // policy.SR.*:
        // policy.*:
        if (not args.skip.policy) {
            if (not args.skip.short_run) {

                state_cache.clear();
                for (eris_time_t t = 0; t < args.data_periods; t++) state_cache.push_back(storage[policy_begins + t]);

                for (auto &d : data) {
                    if (d.applies_to.policy) {
                        if (args.human_readable) output << std::setw(readable_name_width+1) << "policy.SR." + d.name + ":" << " ";
                        else output << ",";

                        output << double_str(d.calculate(storage, policy_begins, policy_begins + args.data_periods - 1), args.double_precision);

                        if (args.human_readable) output << "\n";
                    }
                }
            }

            state_cache.clear();
            for (eris_time_t t = post_policy - args.data_periods; t < post_policy; t++) state_cache.push_back(storage[t]);

            for (auto &d : data) {
                if (d.applies_to.policy) {
                    if (args.human_readable) output << std::setw(readable_name_width+1) << "policy." + d.name + ":" << " ";
                    else output << ",";

                    output << double_str(d.calculate(storage, post_policy - args.data_periods, post_policy - 1), args.double_precision);

                    if (args.human_readable) output << "\n";
                }
            }
        }

        output_mutex.lock();
        std::cout << output.str() << std::endl;
        output_count++;
        output_mutex.unlock();
        output.str("");

        input_lock.lock();
    }
    input_lock.unlock();

    THREAD_DONE
}

int main(int argc, char *argv[]) {
    cmdargs::Data args;
    args.parse(argc, argv);

    // Simulation parameters:
    auto initial_data = data::initial_data_fields();
    auto data = data::data_fields();

    unsigned int longest_name = 0;
    if (args.human_readable) {
        // Figure out long the longest field is to line things up for the human
        for (const auto &d : initial_data) if (d.name.length() + 6 > longest_name) longest_name = d.name.length() + 6;
        for (const auto &d : data) if (d.name.length() + 7 > longest_name) longest_name = d.name.length() + 7;
    }
    else if (not args.no_csv_header) {
        // Write CSV header
        std::ostringstream output;
        output << "source";
        for (const auto &d : initial_data) output << ",param." << d.name;
        for (const auto &d : data) if (d.applies_to.pre) output << ",pre." << d.name;
        if (not args.skip.piracy) {
            if (not args.skip.short_run) for (const auto &d : data) if (d.applies_to.piracy) output << ",piracy.SR." << d.name;
            for (const auto &d : data) if (d.applies_to.piracy) output << ",piracy." << d.name;
        }
        if (not args.skip.policy) {
            if (not args.skip.short_run) for (const auto &d : data) if (d.applies_to.policy) output << ",policy.SR." << d.name;
            for (const auto &d : data) if (d.applies_to.policy) output << ",policy." << d.name;
        }
        output << "\n";
        std::cout << output.str();
    }

    input_count = args.input.size();
    input_it = args.input.cbegin();
    input_it_end = args.input.cend();
    std::vector<std::thread> threads;
    if (args.preload > 0) {
        // Start the preload thread, which loads files into memory:
        THREAD_START
        threads.emplace_back(thr_preload, args);
    }

    if (args.threads == 0) {
        // Not doing threaded loading (aside from, possibly, preloading while parsing):
        THREAD_START // Do this anyway (so the "threads_running" value stays correct)
        thr_parse_file(args, initial_data, data, longest_name);
    }
    else {
        Eigen::initParallel();
        for (unsigned t = 0; t < args.threads; t++) {
            THREAD_START
            threads.emplace_back(thr_parse_file, args, initial_data, data, longest_name);
        }
    }

#ifdef __linux__
    std::unique_lock<std::mutex> thr_lock(thread_mutex);
    while (threads_running > 0) {
        thread_done_cv.wait_for(thr_lock, std::chrono::milliseconds(250), [] { return threads_running == 0; });

        std::lock_guard<std::mutex> lg(output_mutex);
        // Update the process name to something like "crdata [43/123]" (the space and "a" before the [ get
        // eliminated if required--we aren't allowed to set a name longer than 15 characters).
        std::string progress = "[" + std::to_string(processing_counter) + "/" + std::to_string(input_count) + "]";
        std::string name =
            std::string("crdata ").substr(0, std::max<int>(5, 15 - (int) progress.size()))
            + progress;

        prctl(PR_SET_NAME, name.c_str());
    }
#endif

    // Wait for preload and/or parse threads to finish:
    for (auto &th : threads) {
        if (th.joinable()) th.join();
    }

    if (output_count == 0 and not args.only_csv_header) {
        std::cerr << "Error: no usable data sources!\n";
        return 1;
    }
}

