#include "creativity/Creativity.hpp"
#include "creativity/data/Data.hpp"
#include "creativity/state/Storage.hpp"
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
#include <iomanip>
#include <list>

namespace creativity { namespace state { class State; } }

using namespace creativity;
using namespace creativity::state;
using namespace eris;

// When outputting CSV values, we want exact precision, but don't necessarily need the full
// max_digits10 digits to get there.  For example, 0.1 with max_digits10 becomes 0.10000000000000001
// but 0.1 also converts to the numerically identical value.  0.10000000000000002, on the other
// hand, needs every decimal digit of precision.
//
// This function produces tries first at the requested precision, then the requested precision less
// one, then less two; if the subsequent values are numerically identical to the given double value,
// the shortest is returned.
std::string double_str(double d, unsigned precision) {
    double round_trip;
    for (unsigned prec : {precision-2, precision-1}) {
        std::stringstream ss;
        ss.precision(prec);
        ss << d;
        ss >> round_trip;
        if (round_trip == d) { ss << d; return ss.str(); }
    }
    std::stringstream ss;
    ss.precision(precision);
    ss << d;
    return ss.str();
}

unsigned int output_count;
std::mutex output_mutex; // Guards std::out, output_count
decltype(cmdargs::Data::input.cbegin()) input_it, input_it_end;
std::mutex input_it_mutex;

void thr_parse_file(
        const cmdargs::Data &args,
        const std::vector<data::initial_datum> &initial_data,
        const std::vector<data::datum> &data,
        unsigned readable_name_width) {
    input_it_mutex.lock();
    while (input_it != input_it_end) {
        std::string source(*input_it++);
        input_it_mutex.unlock();
        output_mutex.lock();
        std::cerr << "Processing " << source << "\n";
        output_mutex.unlock();

        std::ostringstream output;
        output.precision(args.double_precision);
        auto creativity = Creativity::create();
        // Filename input
        try {
            creativity->fileRead(source);
        }
        catch (std::ios_base::failure&) {
            std::cerr << "Unable to read `" << source << "': " << std::strerror(errno) << "\n";
            continue;
        }
        catch (std::exception &e) {
            std::cerr << "Unable to read `" << source << "': " << e.what() << "\n";
            continue;
        }

        auto locked_storage = creativity->storage();
        auto &storage = *locked_storage.first;

        // Check that the data source has enough data:

        const eris_time_t &piracy_begins = creativity->parameters.piracy_begins,
              &public_sharing_begins = creativity->parameters.public_sharing_begins;
        eris_time_t post_piracy = 0, post_public = 0, post_pre = 0;
        if (not args.skip.piracy)
            post_piracy = (args.skip.public_sharing ? storage.size() : public_sharing_begins);
        if (not args.skip.public_sharing)
            post_public = storage.size();
        post_pre = args.skip.piracy
                ? args.skip.public_sharing
                    ? storage.size()
                    : public_sharing_begins
                : piracy_begins;

#define SKIP_IF(CONDITION, REASON) if (CONDITION) { std::cerr << "Skipping `" << source << "': " << REASON << "\n"; continue; }

        SKIP_IF(args.verify.periods > 0 and (size_t) args.verify.periods != storage.size()-1,
                "simulation periods " << storage.size()-1 << " != " << args.verify.periods);
        SKIP_IF(args.verify.piracy_begins > 0 and (size_t) args.verify.piracy_begins != piracy_begins,
                "simulation piracy_begins " << piracy_begins << " != " << args.verify.piracy_begins);
        SKIP_IF(args.verify.public_sharing_begins > 0 and (size_t) args.verify.public_sharing_begins != public_sharing_begins,
                "simulation public_sharing_begins " << public_sharing_begins << " != " << args.verify.public_sharing_begins);

        SKIP_IF(post_pre-1 <= args.data_periods,
                "simulation \"pre\" periods " << post_pre-2 << " < " << args.data_periods);

        // If we're including both short-run and long-run effects, we need twice as many periods to
        // avoid overlap.
        unsigned periods_needed = args.data_periods;
        if (not args.skip.short_run) periods_needed *= 2;

        if (not args.skip.piracy) {
            SKIP_IF(post_piracy - piracy_begins < periods_needed,
                    "simulation doesn't have enough piracy periods (t=" << piracy_begins << " through t=" << post_piracy-1 << ")");
        }
        if (not args.skip.public_sharing) {
            SKIP_IF(post_public - public_sharing_begins < periods_needed,
                    "simulation doesn't have enough public sharing periods (t=" << public_sharing_begins << " through t=" << post_public-1 << ")");
        }
#undef SKIP_IF

        if (args.human_readable) output << "\n\n" << source << "\n==========\n";
        else output << data::csv_fix(source);
        for (auto &d : initial_data) {
            if (args.human_readable) output << std::setw(readable_name_width+1) << d.name + ":" << " ";
            else output << ",";

            if (d.calc_double) output << double_str(d.calc_double(creativity->parameters), args.double_precision);
            else output << d.calc_int(creativity->parameters);

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

        // public.SR.*:
        // public.*:
        if (not args.skip.public_sharing) {
            if (not args.skip.short_run) {

                state_cache.clear();
                for (eris_time_t t = 0; t < args.data_periods; t++) state_cache.push_back(storage[public_sharing_begins + t]);

                for (auto &d : data) {
                    if (d.applies_to.public_sharing) {
                        if (args.human_readable) output << std::setw(readable_name_width+1) << "public.SR." + d.name + ":" << " ";
                        else output << ",";

                        output << double_str(d.calculate(storage, public_sharing_begins, public_sharing_begins + args.data_periods - 1), args.double_precision);

                        if (args.human_readable) output << "\n";
                    }
                }
            }

            state_cache.clear();
            for (eris_time_t t = post_public - args.data_periods; t < post_public; t++) state_cache.push_back(storage[t]);

            for (auto &d : data) {
                if (d.applies_to.public_sharing) {
                    if (args.human_readable) output << std::setw(readable_name_width+1) << "public." + d.name + ":" << " ";
                    else output << ",";

                    output << double_str(d.calculate(storage, post_public - args.data_periods, post_public - 1), args.double_precision);

                    if (args.human_readable) output << "\n";
                }
            }
        }

        output_mutex.lock();
        std::cout << output.str() << std::endl;
        output_count++;
        output_mutex.unlock();
        output.str("");

        input_it_mutex.lock();
    }
    input_it_mutex.unlock();
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
        if (not args.skip.public_sharing) {
            if (not args.skip.short_run) for (const auto &d : data) if (d.applies_to.public_sharing) output << ",public.SR." << d.name;
            for (const auto &d : data) if (d.applies_to.public_sharing) output << ",public." << d.name;
        }
        output << "\n";
        std::cout << output.str();
    }

    input_it = args.input.cbegin();
    input_it_end = args.input.cend();
    if (args.threads == 0) {
        thr_parse_file(args, initial_data, data, longest_name);
    }
    else {
        std::vector<std::thread> threads;
        for (unsigned t = 0; t < args.threads; t++) {
            threads.emplace_back(thr_parse_file, args, initial_data, data, longest_name);
        }
        for (auto &th : threads) {
            if (th.joinable()) th.join();
        }
    }

    if (output_count == 0 and not args.only_csv_header) {
        std::cerr << "Error: no usable data sources!\n";
        return 1;
    }
}

