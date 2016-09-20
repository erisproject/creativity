#include <eris/random/rng.hpp>
#include <eris/random/util.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <eris/random/truncated_normal_distribution.hpp>
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

using namespace eris;

const std::string help_message = u8R"(
This script is a wrapper intended for use with creativity-cli or creativity-gui (but usable with any
executable) that allows argument values to be drawn from distributions.  In particular, it looks for
any arguments of the following forms:

    U[a,b] - draws a uniformly distributed double from [a, b)
    iU[a,b] - draws a uniformly distributed integer from [a,b]
    N(m,sd) - draws a normally distributed double with mean m and standard deviation sd
    N(m,sd)[a,b] - draws a normally distributed double, truncated to the range [a, b].
                   One of 'a' or 'b' can be omitted to specify a one-sided truncation.
    N(m,sd)+ - equivalent to N(m,sd)[0,].
    iN(m,sd) - like N(m,sd), but rounds the drawn value to the nearest integer.
    iN(m,sd)[a,b] - Truncated version of iN(m,sd).
    {a,b,c,...} - selects one of the comma-separated strings "a", "b", "c", etc. with equal probability.

Any other argument is passed through as is.  At least one argument matching the above must be included.

Example:

    ./creativity-random ./creativity-cli -D 2 -r 'iN(100,25)[50,150]' -f 'U[0.01,0.25]' -Q 'N(0,1)+'

which would take the appropriate draws and could, for example, execute the following:

    ./creativity-cli -D 2 -r 105 -f 0.19650807491555314 -Q 0.86839032695624874

Note that if ERIS_RNG_SEED is set, the given seed is used (independently) for *both* the
creativity-random and creativity-cli processes.
)";


const std::string
    re_double("[-+]?[0-9]*\\.?[0-9]+(?:e[-+]?[0-9]+)?"),
    re_double_pos("\\+?[0-9]*\\.?[0-9]+(?:e[-+]?[0-9]+)?)"),
    re_double_inf("[-+]?(?:[0-9]*\\.?[0-9]+(?:e[-+]?[0-9]+)?|inf(?:inity)?)"),
    re_int("[-+]?[0-9]+"),
    re_int_pos("\\+?[0-9]+");

std::vector<std::pair<std::regex, std::function<std::string(const std::smatch&)>>> load_patterns() {
    decltype(load_patterns()) callbacks;
    // U[a,b]
    callbacks.emplace_back(std::regex("^U[\\[(](" + re_double + ")\\s*,\\s*(" + re_double + ")[\\])]$", std::regex::icase),
        [](const std::smatch &m) -> std::string {
            double a = std::stod(m[1]);
            double b = std::stod(m[2]);
            if (b <= a) throw std::logic_error("Error in input: `" + m[0].str() + "' is invalid (max <= min)");
            std::ostringstream result;
            result << std::setprecision(std::numeric_limits<double>::max_digits10)
                << random::runiform(a, b);
            return result.str();
        });
    // iU[a,b]
    callbacks.emplace_back(std::regex("^iU[\\[(](" + re_int + ")\\s*,\\s*(" + re_int + ")[\\])]$", std::regex::icase),
        [](const std::smatch &m) -> std::string {
            long a = std::stol(m[1]);
            long b = std::stol(m[2]);
            if (b <= a) throw std::logic_error("Error in input: `" + m[0].str() + "' is invalid (max <= min)");
            return std::to_string(boost::random::uniform_int_distribution<long>(a, b)(random::rng()));
        });
    // All the normal variants
    callbacks.emplace_back(std::regex("^(i)?N\\((" + re_double + ")\\s*,\\s*(" + re_double + ")\\)" + 
                "(?:(\\+)|\\[(" + re_double_inf + ")?\\s*,\\s*(" + re_double_inf + ")?\\])?$", std::regex::icase),
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
            if (sd < 0) throw std::logic_error("Error in input: `" + m[0].str() + "' has invalid negative standard deviation");
            if (ubound < lbound) throw std::logic_error("Error in input: `" + m[0].str() + "' has impossible-to-satisfy constraint");
            std::function<double()> draw_norm;
            if (std::isfinite(ubound) or std::isfinite(lbound)) {
                draw_norm = [&] () -> double { return random::truncated_normal_distribution<double>(mean, sd, lbound, ubound)(random::rng()); };
            }
            else {
                draw_norm = [&] () -> double { return random::normal_distribution<double>(mean, sd)(random::rng()); };
            }
            double draw;
            int draws = 0;
            do {
                if (draws++ > 1000) throw std::runtime_error("Draw failure: `" + m[0].str() + "' constraints were not satisfied within 1000 draws");
                draw = draw_norm();
                if (draw < 0 and mean == 0 and lbound == 0) draw = -draw; // Optimize special case (absolute value of 0-mean normal)
                if (round) draw = std::round(draw);
            } while (draw < lbound or draw > ubound);

            if (round) return std::to_string(std::lround(draw));

            std::ostringstream result;
            result << std::setprecision(std::numeric_limits<double>::max_digits10) << draw;
            return result.str();
        });

    // Random string from comma-separated strings
    callbacks.emplace_back(std::regex("\\{(.*)\\}"), [](const std::smatch &m) -> std::string {
        std::string result;
        std::istringstream iss;
        iss.str(m[1].str());
        std::string item;
        uint_fast32_t count = 0;
        while (std::getline(iss, item, ',')) {
            if (eris::random::rcoin(1.0 / ++count))
                result = item;
        }

        return result;
    });

    return callbacks;
}

int main (int argc, const char* argv[]) {
    // This vector stores the patterns and associated lambdas to call, and will be tried in order
    // until a match is found.
    auto try_match = load_patterns();

    // The first argument is the program to execute, but make sure it is specified and isn't set to --help
    bool help = (argc < 3 or std::string(argv[1]) == "--help"); // ... unless the first argument is --help
    bool found = false;
    std::vector<std::string> args;
    if (argc >= 3) args.reserve(argc-2);
    for (int i = 2; i < argc; i++) args.push_back(argv[i]);
    for (auto &arg : args) {
        if (arg == "--help") help = true;
        else for (auto &p : try_match) {
            std::smatch match_res;
            if (std::regex_match(arg, match_res, p.first)) {
                found = true;
                arg = p.second(match_res);
                break;
            }
        }
    }

    // --help gets passed through, but we *also* handle it by printing help for creativity-random
    if (help or not found) {
        std::cout << "USAGE: " << argv[0] << " PROGRAM ARG [ARG ...]\n" << help_message << "\n";
        if (not found) {
            std::cerr << "Aborting because no random signatures found in argument list!\n";
            exit(1);
        }
    }

    const char **exec_argv = new const char*[2 + args.size()]; // +1 for the executable, +1 for the terminating nullptr
    {
        int i = 0;
        std::cout << "Executing " << argv[1];
        exec_argv[i++] = argv[1];
        for (auto &arg : args) {
            exec_argv[i++] = arg.c_str();
            std::cout << " " << arg;
        }
        exec_argv[i] = nullptr;
        std::cout << "\n";
    }

    execvp(exec_argv[0], const_cast<char**>(exec_argv));

    // Getting here means exec failed!
    std::cerr << "Failed to execute `" << exec_argv[0] << "': " << std::strerror(errno) << "\n";

    std::exit(errno);
}
