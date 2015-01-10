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
    unsigned int pre_piracy_periods = 50,
                 new_piracy_periods = 50,
                 post_piracy_periods = 50,
                 double_precision = std::numeric_limits<double>::digits10;

    std::vector<std::string> input;
};
args parseArgs(int argc, const char *argv[]) {
    TCLAP::CmdLine cmd("Eris-based creativity simulation model - data generator", ' ', "0.0.1");

    args a;

    TCLAP::ValueArg<long> piracy_begins_arg("P", "verify-piracy-begins", "Only use given files with a piracy_begins value matching the given argument", false, a.check_piracy_begins, "integer");
    TCLAP::ValueArg<long> periods_arg("T", "verify-periods", "Only use given files with the last time period matching the given argument", false, a.check_periods, "integer");
    TCLAP::ValueArg<unsigned int> pre_piracy_periods_arg("t", "pre", "Use the last 'p' pre-piracy periods to calculate pre-piracy values", true, a.pre_piracy_periods, "integer");
    TCLAP::ValueArg<unsigned int> new_piracy_periods_arg("p", "new", "Use the first 'n' post-piracy periods to calculate the new-piracy values", true, a.new_piracy_periods, "integer");
    TCLAP::ValueArg<unsigned int> post_piracy_periods_arg("z", "post", "Use the last 't' periods to calculate the post-piracy values", true, a.post_piracy_periods, "integer");
    TCLAP::ValueArg<unsigned int> double_precision_arg("", "precision", "Specifies the precision level for double values. The default, " + std::to_string(a.double_precision) + ", is the minimum required to exactly represent all possible double values", false, a.double_precision, "integer");

    cmd.add(double_precision_arg);
    cmd.add(piracy_begins_arg);
    cmd.add(periods_arg);
    cmd.add(post_piracy_periods_arg);
    cmd.add(new_piracy_periods_arg);
    cmd.add(pre_piracy_periods_arg);

    //TCLAP::UnlabeledMultiArg<std::string> input("CRSTATE files", "One or more .crstate files to generate statistics from", true, "filenames");
    //cmd.add(input);
    TCLAP::UnlabeledMultiArg<std::string> input("input-files", "One or more .crstate files to generate statistics from", true, "DATASOURCE");
    cmd.add(input);

    cmd.parse(argc, argv);

    a.check_piracy_begins = piracy_begins_arg.getValue();
    a.check_periods = periods_arg.getValue();
    a.pre_piracy_periods = pre_piracy_periods_arg.getValue();
    a.new_piracy_periods = new_piracy_periods_arg.getValue();
    a.post_piracy_periods = post_piracy_periods_arg.getValue();
    a.double_precision = double_precision_arg.getValue();
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
    output << "source";
    for (const auto &d : initial_data) output << "," << d.name;
    for (const auto &d : data) if (!d.piracy_only) output << "," << "pre_" << d.name;
    for (const auto &d : data) output << "," << "new_" << d.name;
    for (const auto &d : data) output << "," << "post_" << d.name;
    output << "\n";

    for (const auto &source : a.input) {
        auto creativity = Creativity::create();
        if (source.substr(0, 13) == "postgresql://" or source.substr(0, 11) == "postgres://") {
            try {
                creativity->pgsql(source, true /*read-only*/);
            }
            catch (std::exception &e) {
                std::cerr << "Unable to connect to database `" << source << "': " << e.what() << "\n\n";
                continue;
            }
        }
        else {
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
        }

        auto locked_storage = creativity->storage();
        auto &storage = *locked_storage.first;

        // Check that the data source has enough data:

#define SKIP_IF(CONDITION, REASON) if (CONDITION) { std::cerr << "Skipping `" << source << "': " << REASON << "\n"; continue; }

        SKIP_IF(a.check_periods > 0 and (size_t) a.check_periods != storage.size()-1,
                "simulation periods " << storage.size()-1 << " != " << a.check_periods);
        SKIP_IF(a.check_piracy_begins >= 0 and (size_t) a.check_piracy_begins != creativity->parameters.piracy_begins,
                "simulation piracy_begins " << creativity->parameters.piracy_begins << " != " << a.check_piracy_begins);
        SKIP_IF(creativity->parameters.piracy_begins <= a.pre_piracy_periods,
                "simulation per-piracy periods " << creativity->parameters.piracy_begins << " != " << a.pre_piracy_periods);
        SKIP_IF(storage.size() <= creativity->parameters.piracy_begins,
                "simulation doesn't have any piracy periods (T=" << storage.size()-1 << ", P=" << creativity->parameters.piracy_begins << ")");
        // e.g. size = 200 means we have T=0 through T=199; if P=150, then there are 200-150
        // (==199-150+1) = 50 piracy periods (150 through 199), so we could allow any new+post <=
        // 50.  Or equivalently, size must be at least piracy_begins + new + post:
        SKIP_IF(storage.size() < creativity->parameters.piracy_begins + a.new_piracy_periods + a.post_piracy_periods,
                "simulation has insufficient post-piracy periods (need " << a.new_piracy_periods << " (new) + " << a.post_piracy_periods <<
                " (post) but only have " << storage.size() - creativity->parameters.piracy_begins << " post-piracy periods");

#undef SKIP_IF

        const eris_time_t &piracy_begins = creativity->parameters.piracy_begins;

        output << data::csv_escape(source);
        for (auto &d : initial_data) {
            if (d.calc_double) output << "," << d.calc_double(creativity->parameters);
            else output << "," << d.calc_int(creativity->parameters);
        }

        // Keep a reference to all of the periods, so that the following don't end up continually
        // reloading them from disk/database.
        std::list<std::shared_ptr<const State>> state_cache;
        for (eris_id_t t = piracy_begins - a.pre_piracy_periods; t < piracy_begins + a.new_piracy_periods; t++)
            state_cache.push_back(storage[t]);
        for (eris_id_t t = storage.size() - a.post_piracy_periods; t < storage.size(); t++)
            state_cache.push_back(storage[t]);

        for (auto &d : data) {
            if (not d.piracy_only)
                output << "," << d.calculate(storage, piracy_begins - a.pre_piracy_periods, piracy_begins - 1);
        }

        for (auto &d : data)
            output << "," << d.calculate(storage, piracy_begins, piracy_begins + a.new_piracy_periods - 1);

        for (auto &d : data)
            output << "," << d.calculate(storage, storage.size() - a.post_piracy_periods, storage.size() - 1);

        output << "\n";
        output_count++;
    }

    if (output_count > 0)
        std::cout << output.str();
    else
        std::cerr << "Error: no usable data sources!\n";
}

