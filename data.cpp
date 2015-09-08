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

int main(int argc, char *argv[]) {
    cmdargs::Data args;
    args.parse(argc, argv);

    // Simulation parameters:
    auto initial_data = data::initial_data_fields();
    auto data = data::data_fields();

    std::ostringstream output;
    int output_count = 0;
    output << std::setprecision(args.double_precision);
    unsigned int longest_name = 0;
    if (not args.human_readable) {
        output << "source";
        for (const auto &d : initial_data) output << ",param_" << d.name;
        for (const auto &d : data) if (d.applies_to.pre) output << ",pre_" << d.name;
        for (const auto &d : data) if (d.applies_to.piracy) output << ",piracy_" << d.name;
        for (const auto &d : data) if (d.applies_to.public_sharing) output << "," << "public_" << d.name;
        output << "\n";
        std::cout << output.str();
        output.str("");
    }
    else {
        for (const auto &d : initial_data) if (d.name.length() + 6 > longest_name) longest_name = d.name.length() + 6;
        for (const auto &d : data) if (d.name.length() + 7 > longest_name) longest_name = d.name.length() + 7;
    }

    for (const auto &source : args.input) {
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
        if (not args.skip.piracy) {
            SKIP_IF(post_piracy - piracy_begins < args.data_periods,
                    "simulation doesn't have enough piracy periods (t=" << piracy_begins << " through t=" << post_piracy-1 << ")");
        }
        if (not args.skip.public_sharing) {
            SKIP_IF(post_public - public_sharing_begins < args.data_periods,
                    "simulation doesn't have enough public sharing periods (t=" << public_sharing_begins << " through t=" << post_public-1 << ")");
        }
#undef SKIP_IF

        if (args.human_readable) output << "\n\n" << source << "\n==========\n";
        else output << data::csv_escape(source);
        for (auto &d : initial_data) {
            if (args.human_readable) output << std::setw(longest_name+1) << d.name + ":" << " ";
            else output << ",";

            if (d.calc_double) output << d.calc_double(creativity->parameters);
            else output << d.calc_int(creativity->parameters);

            if (args.human_readable) output << "\n";
        }

        // Keep a reference to all of the periods, so that the following don't end up continually
        // reloading them from disk/database.
        std::list<std::shared_ptr<const State>> state_cache;

        for (eris_time_t t = post_pre    - args.data_periods; t < post_pre;    t++) state_cache.push_back(storage[t]);
        for (eris_time_t t = post_piracy - args.data_periods; t < post_piracy; t++) state_cache.push_back(storage[t]);
        for (eris_time_t t = post_public - args.data_periods; t < post_public; t++) state_cache.push_back(storage[t]);

        // pre_*:
        for (auto &d : data) {
            if (d.applies_to.pre) {
                if (args.human_readable) output << std::setw(longest_name+1) << "pre_" + d.name + ":" << " ";
                else output << ",";

                output << d.calculate(storage, post_pre - args.data_periods, post_pre - 1);

                if (args.human_readable) output << "\n";
            }
        }

        // piracy_*:
        if (not args.skip.piracy) {
            for (auto &d : data) {
                if (d.applies_to.piracy) {
                    if (args.human_readable) output << std::setw(longest_name+1) << "piracy_" + d.name + ":" << " ";
                    else output << ",";

                    output << d.calculate(storage, post_piracy - args.data_periods, post_piracy - 1);

                    if (args.human_readable) output << "\n";
                }
            }
        }

        // public_*:
        if (not args.skip.public_sharing) {
            for (auto &d : data) {
                if (d.applies_to.public_sharing) {
                    if (args.human_readable) output << std::setw(longest_name+1) << "public_" + d.name + ":" << " ";
                    else output << ",";

                    output << d.calculate(storage, post_public - args.data_periods, post_public - 1);

                    if (args.human_readable) output << "\n";
                }
            }
        }

        std::cout << output.str() << std::endl;
        output.str("");
        output_count++;
    }

    if (output_count == 0)
        std::cerr << "Error: no usable data sources!\n";
}

