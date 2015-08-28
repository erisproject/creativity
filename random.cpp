#include <eris/Random.hpp>
#include <regex>
#include <functional>
#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>
#include <cmath>
#include <cstring>
#include <cerrno>
#include <cstdlib>
extern "C" {
#include <unistd.h>
}

const std::string help_message = u8R"(
This script is a wrapper around creativity-cli that allows argument values to be drawn from
distributions.  In particular, it looks for any arguments of the following forms:

    U[a,b] - draws a uniformly distributed double from [a, b)
    iU[a,b] - draws a uniformly distributed integer from [a,b]
    N(m,sd) - draws a normally distributed double with mean m and standard deviation sd
    N(m,sd)[a,b] - draws a normally distributed double, truncated to the range [a, b].  The value is
                   redrawn until a value in the range is found.  'a' or 'b' can be omitted to not
                   specify an lower or upper bound, respectively.
    N(m,sd)+ - equivalent to N(m,sd)[0,].
    iN(m,sd) - like N(m,sd), but rounds the drawn value to the nearest integer.
    iN(m,sd)[a,b] - Like N(m,sd), but rounds the drawn value to the nearest integer, and repeats
                    until a (rounded) value in [a,b] is found.

Any other argument is passed through as is.

Example:

    ./creativity-random -D 2 -r 'iN(100,25)[50,150]' -f 'U[0.01,0.25]' -Q 'N(0,1)+'

which would take the appropriate draws and could, for example, execute the following:

    ./creativity-cli -D 2 -r 105 -f 0.19650807491555314 -Q 0.86839032695624874

Note that if ERIS_RNG_SEED is set, the given seed is used (independently) for *both* the
creativity-random and creativity-cli processes.
)";


const std::string
    re_double("[-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?"),
    re_double_pos("\\+?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?)"),
    re_double_inf("[-+]?(?:[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?|[iI][nN][fF](?:[iI][nN][iI][tT][yY])?)"),
    re_int("[-+]?[0-9]+"),
    re_int_pos("\\+?[0-9]+");

std::vector<std::pair<std::regex, std::function<std::string(const std::smatch&)>>> load_patterns() {
    decltype(load_patterns()) callbacks;
    // U[a,b]
    callbacks.emplace_back(std::regex("^U[\\[(](" + re_double + ")\\s*,\\s*(" + re_double + ")[\\])]$"),
        [](const std::smatch &m) -> std::string {
            double a = std::stod(m[1]);
            double b = std::stod(m[2]);
            if (b <= a) throw std::logic_error("Error in input: `" + m[0].str() + "' is invalid (max <= min)");
            std::ostringstream result;
            result << std::setprecision(std::numeric_limits<double>::max_digits10)
                << std::uniform_real_distribution<double>(a, b)(eris::Random::rng());
            return result.str();
        });
    // iU[a,b]
    callbacks.emplace_back(std::regex("^iU[\\[(](" + re_int + ")\\s*,\\s*(" + re_int + ")[\\])]$"),
        [](const std::smatch &m) -> std::string {
            long a = std::stol(m[1]);
            long b = std::stol(m[2]);
            if (b <= a) throw std::logic_error("Error in input: `" + m[0].str() + "' is invalid (max <= min)");
            return std::to_string(std::uniform_int_distribution<long>(a, b)(eris::Random::rng()));
        });
    // All the normal variants
    callbacks.emplace_back(std::regex("^(i)?N\\((" + re_double + ")\\s*,\\s*(" + re_double + ")\\)" + 
                "(?:(\\+)|\\[(" + re_double_inf + ")?\\s*,\\s*(" + re_double_inf + ")?\\])?$"),
        [](const std::smatch &m) -> std::string {
            bool round = m[1].matched;
            double mean = std::stod(m[2]);
            double sd = std::stod(m[3]);
            double lbound, ubound;
            if (m[4].matched) {
                lbound = 0;
                ubound = std::numeric_limits<double>::infinity();
            }
            else {
                lbound = m[5].matched ? std::stod(m[5]) : -std::numeric_limits<double>::infinity();
                ubound = m[6].matched ? std::stod(m[6]) : std::numeric_limits<double>::infinity();
            }
            if (ubound <= lbound) throw std::logic_error("Error in input: `" + m[0].str() + "' has impossible-to-satisfy constraint");
            if (sd < 0) throw std::logic_error("Error in input: `" + m[0].str() + "' has invalid negative standard deviation");
            std::normal_distribution<double> rnorm(mean, sd);
            double draw;
            int draws = 0;
            do {
                if (draws++ > 1000) throw std::runtime_error("Draw failure: `" + m[0].str() + "' constraints were not satisfied within 1000 draws");
                draw = rnorm(eris::Random::rng());
                if (draw < 0 and mean == 0 and lbound == 0) draw = -draw; // Optimize special case (absolute value of 0-mean normal)
                if (round) draw = std::round(draw);
            } while (draw < lbound or draw > ubound);

            if (round) return std::to_string(std::lround(draw));

            std::ostringstream result;
            result << std::setprecision(std::numeric_limits<double>::max_digits10) << draw;
            return result.str();
        });
    
    return callbacks;
}

int main (int argc, char* argv[]) {
    // Build the path to creativity-cli by replacing "-random" with "-cli" in argv[0].
    // Note that there might be a suffix (e.g. .exe for Windows users)
    std::string cli = std::regex_replace(argv[0], std::regex("-random(?!.*-random)"), "-cli");

    // This vector stores the patterns and associated lambdas to call, and will be tried in order
    // until a match is found.

    auto try_match = load_patterns();

    bool help = false;
    bool found = false;
    std::vector<std::string> args;
    args.reserve(argc-1);
    for (int i = 1; i < argc; i++) args.push_back(argv[i]);
    for (auto &arg : args) {
        if (arg == "--help") help = true;
        for (auto &p : try_match) {
            std::smatch match_res;
            if (std::regex_match(arg, match_res, p.first)) {
                found = true;
                arg = p.second(match_res);
                break;
            }
        }
    }

    std::vector<char*> cli_argv;
    cli_argv.reserve(2 + args.size());
    cli_argv.push_back(const_cast<char*>(cli.c_str()));
    for (auto &arg : args) cli_argv.push_back(const_cast<char*>(arg.c_str()));
    cli_argv.push_back(nullptr);

    // --help gets passed through, but we *also* handle it by printing help for creativity-random
    if (help or not found) {
        std::cout << "USAGE: " << argv[0] << " ARG ...\n" << help_message << "\n";
        if (not found) {
            std::cout << "Aborting because no random signatures found in argument list!\n";
            exit(1);
        }
    }

    std::cout << "Executing";
    for (auto &argv : cli_argv) if (argv) std::cout << " " << argv;
    std::cout << "\n";

    execvp(cli_argv[0], cli_argv.data());

    // Getting here means exec failed!
    std::cerr << "Failed to execute: " << std::strerror(errno) << "\n";

    std::exit(errno);
}
