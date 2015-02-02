#include "creativity/Creativity.hpp"
#include "creativity/data/Data.hpp"
#include <tclap/CmdLine.h>
#include <cstdio>
#include <algorithm>

using namespace creativity;
using namespace creativity::state;
using namespace eris;

struct args {
    long check_piracy_begins = -1,
         check_periods = -1;
    unsigned int periods = 25,
                 skip_piracy = 1,
                 double_precision = std::numeric_limits<double>::digits10;
    bool human_readable = false;

    std::vector<std::string> input;
};
args parseArgs(int argc, const char *argv[]) {
    TCLAP::CmdLine cmd("Eris-based creativity simulation model - data generator", ' ', "0.0.1");

    args a;

    TCLAP::ValueArg<long> verify_piracy_begins_arg("P", "verify-piracy-begins", "Only use given files with a piracy_begins value matching the given argument", false, a.check_piracy_begins, "integer");
    TCLAP::ValueArg<long> verify_periods_arg("T", "verify-periods", "Only use given files with the last time period matching the given argument", false, a.check_periods, "integer");
    TCLAP::ValueArg<unsigned int> periods_arg("t", "periods", "Specifies the number of periods to use for calculating pre, new, and post-piracy data", true, a.periods, "integer");
    TCLAP::ValueArg<unsigned int> skip_piracy_arg("s", "skip-piracy", "Skips the first 's' piracy periods.  Defaults to 1, as the first piracy period often involves many mispredicted quality draws.", false, a.skip_piracy, "integer");
    TCLAP::ValueArg<unsigned int> double_precision_arg("", "precision", "Specifies the precision level for double values. The default, " + std::to_string(a.double_precision) + ", is the minimum required to exactly represent all possible double values", false, a.double_precision, "integer");
    TCLAP::SwitchArg human_readable("H", "human-readable", "Produce output in human-readable format.  The default (when this is not specified) outputs comma-separated values.  This also changes the default value for --precision to 8.");

    cmd.add(human_readable);
    cmd.add(double_precision_arg);
    cmd.add(verify_piracy_begins_arg);
    cmd.add(verify_periods_arg);
    cmd.add(skip_piracy_arg);
    cmd.add(periods_arg);

    TCLAP::UnlabeledMultiArg<std::string> input("input-files", "One or more .crstate files to generate statistics from", true, "filenames");
    cmd.add(input);

    cmd.parse(argc, argv);

    a.check_piracy_begins = verify_piracy_begins_arg.getValue();
    a.check_periods = verify_periods_arg.getValue();
    a.periods = periods_arg.getValue();
    a.skip_piracy = skip_piracy_arg.getValue();
    a.human_readable = human_readable.getValue();
    a.double_precision = (double_precision_arg.isSet() or not a.human_readable) ? double_precision_arg.getValue() : 8;
    a.input = input.getValue();

    return a;
}


int main(int argc, const char *argv[]) {
    args a = parseArgs(argc, argv);

    // Simulation parameters:
    auto initial_data = data::initial_data_fields();
    auto data = data::data_fields();

    std::ostringstream output;
    int output_count = 0;
    output << std::setprecision(a.double_precision);
    unsigned int longest_name = 0;
    if (not a.human_readable) {
        output << "source";
        for (const auto &d : initial_data) output << "," << d.name;
        for (const auto &d : data) if (!d.piracy_only) output << "," << "pre_" << d.name;
        for (const auto &d : data) output << "," << "new_" << d.name;
        for (const auto &d : data) output << "," << "post_" << d.name;
        output << "\n";
        std::cout << output.str();
        output.str("");
    }
    else {
        for (const auto &d : initial_data) if (d.name.length() > longest_name) longest_name = d.name.length();
        for (const auto &d : data) if (d.name.length() + 5 > longest_name) longest_name = d.name.length() + 5;
    }

    for (const auto &source : a.input) {
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

#define SKIP_IF(CONDITION, REASON) if (CONDITION) { std::cerr << "Skipping `" << source << "': " << REASON << "\n"; continue; }

        SKIP_IF(a.check_periods > 0 and (size_t) a.check_periods != storage.size()-1,
                "simulation periods " << storage.size()-1 << " != " << a.check_periods);
        SKIP_IF(a.check_piracy_begins >= 0 and (size_t) a.check_piracy_begins != creativity->parameters.piracy_begins,
                "simulation piracy_begins " << creativity->parameters.piracy_begins << " != " << a.check_piracy_begins);
        SKIP_IF(creativity->parameters.piracy_begins <= a.periods,
                "simulation per-piracy periods " << creativity->parameters.piracy_begins << " != " << a.periods);
        SKIP_IF(storage.size() <= creativity->parameters.piracy_begins,
                "simulation doesn't have any piracy periods (T=" << storage.size()-1 << ", P=" << creativity->parameters.piracy_begins << ")");
        // e.g. size = 200 means we have T=0 through T=199; if P=150, then there are 200-150
        // (==199-150+1) = 50 piracy periods (150 through 199), so we could allow any new+post <=
        // 50.  Or equivalently, size must be at least piracy_begins + new + post:
        SKIP_IF(storage.size() < creativity->parameters.piracy_begins + a.skip_piracy + 2*a.periods,
                "simulation has insufficient post-piracy periods: need " << a.skip_piracy + 2*a.periods << " ["<< a.skip_piracy << " skip, " << a.periods << " new, " << a.periods <<
                " post], but only have " << storage.size() - creativity->parameters.piracy_begins << " post-piracy periods");

#undef SKIP_IF

        const eris_time_t &piracy_begins = creativity->parameters.piracy_begins;

        if (a.human_readable) output << "\n\n" << source << "\n==========\n";
        else output << data::csv_escape(source);
        for (auto &d : initial_data) {
            if (a.human_readable) output << std::setw(longest_name+1) << d.name + ":" << " ";
            else output << ",";

            if (d.calc_double) output << d.calc_double(creativity->parameters);
            else output << d.calc_int(creativity->parameters);

            if (a.human_readable) output << "\n";
        }

        // Keep a reference to all of the periods, so that the following don't end up continually
        // reloading them from disk/database.
        std::list<std::shared_ptr<const State>> state_cache;
        for (eris_id_t t = piracy_begins - a.periods; t < piracy_begins; t++)
            state_cache.push_back(storage[t]);
        for (eris_id_t t = piracy_begins + a.skip_piracy; t < piracy_begins + a.skip_piracy + a.periods; t++)
            state_cache.push_back(storage[t]);
        for (eris_id_t t = storage.size() - a.periods; t < storage.size(); t++)
            state_cache.push_back(storage[t]);

        // pre_*:
        for (auto &d : data) {
            if (not d.piracy_only) {
                if (a.human_readable) output << std::setw(longest_name+1) << "pre_" + d.name + ":" << " ";
                else output << ",";

                output << d.calculate(storage, piracy_begins - a.periods, piracy_begins - 1);

                if (a.human_readable) output << "\n";
            }
        }
        std::cout << "skip: " << a.skip_piracy << "\n";
        std::cout << "skip: " << a.periods << "\n";

        // new_*:
        for (auto &d : data) {
            if (a.human_readable) output << std::setw(longest_name+1) << "new_" + d.name + ":" << " ";
            else output << ",";

            output << d.calculate(storage, piracy_begins + a.skip_piracy, piracy_begins + a.skip_piracy + a.periods - 1);

            if (a.human_readable) output << "\n";
        }

        // post_*:
        for (auto &d : data) {
            if (a.human_readable) output << std::setw(longest_name+1) << "post_" + d.name + ":" << " ";
            else output << ",";

            output << d.calculate(storage, storage.size() - a.periods, storage.size() - 1);

            if (a.human_readable) output << "\n";
        }

        output << "\n";
        std::cout << output.str();
        output.str("");
        output_count++;
    }

    if (output_count == 0)
        std::cerr << "Error: no usable data sources!\n";
}

