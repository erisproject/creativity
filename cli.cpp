#include "creativity/Creativity.hpp"
#include "creativity/state/Storage.hpp"
#ifndef CREATIVITY_SKIP_PGSQL
#include "creativity/state/PsqlStorage.hpp"
#endif
#include <eris/Simulation.hpp>
#include <eris/Random.hpp>
#include <functional>
#include <iostream>
#include <iomanip>
#include <regex>
#include <Eigen/Core>
#include <tclap/CmdLine.h>
#include <sys/stat.h>
#include <typeinfo>

using namespace creativity;
using namespace creativity::state;
using namespace eris;
using namespace Eigen;

// For some unknown reason, TCLAP adds arguments to its internal list in reverse order, so the
// generated --help output will be in reverse order from the order arguments are added in the code.
// This subclass is here to expose the list so that it can be hacked around (i.e. reversed).
//
// Call startReverseHack() just before adding the first argument, then call reverseHack() after
// adding the last argument.
class CmdLineHack : public TCLAP::CmdLine {
    private:
        decltype(_argList.begin()) rev_end_;
    public:
        using CmdLine::CmdLine;
        // Records the current beginning of the list as the first element not to be reversed
        void startReverseHack() { rev_end_ = _argList.begin(); }
        // Reverse the argument list from the beginning up to what the beginning was when
        // startReverseHack() was called.
        void reverseHack() { std::reverse(_argList.begin(), rev_end_); }
};

// Returns a double converted to a string, trimming off insignificant 0s.
template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
std::string output_string(T v) {
    return std::to_string(v);
}

template <>
std::string output_string(double v) {
    std::string z = std::to_string(v);
    std::regex trim_trailing("(\\.\\d*?)0+$"), trim_empty_dot("\\.$");
    return std::regex_replace(std::regex_replace(z, trim_trailing, "$1"), trim_empty_dot, "");
}

/// Constraint that takes a numeric minimum (and optionally a maximum); the argument must be >= min
/// and <= max.
template <class T, class = typename std::enable_if<std::is_arithmetic<T>::value>::type>
class RangeConstraint : public TCLAP::Constraint<T> {
    public:
        RangeConstraint() = default; // Bounded only by data type limits
        RangeConstraint(T minimum, T maximum, bool strict = false) : min(minimum), max(maximum), strict(strict) {}
        static RangeConstraint LE(T maximum) { return RangeConstraint(LOWEST, maximum); }
        static RangeConstraint GE(T minimum) { return RangeConstraint(minimum, BIGGEST); }
        static RangeConstraint LT(T maximum) { return RangeConstraint(LOWEST, maximum, true); }
        static RangeConstraint GT(T minimum) { return RangeConstraint(minimum, BIGGEST, true); }
        virtual std::string description() const override {
            if (max == BIGGEST)
                return (strict ? "Value must be greater than " : "Value must be at least ") + output_string(min);
            else if (min == LOWEST)
                return (strict ? "Value must be less than " : "Value must be at most ") + output_string(min);
            else
                return (strict ? "Value must be strictly between " : "Value must be between ")
                    + output_string(min) + " and " + output_string(max);
        }
        virtual std::string shortID() const override {
            std::string type =
                std::is_floating_point<T>::value ? u8"𝑹" :
                std::is_integral<T>::value ? (std::is_unsigned<T>::value ? u8"𝑵" : u8"𝒁") :
                "v"; // fallback
            bool show_min = (min != LOWEST or LOWEST == 0),
                 show_max = (max != BIGGEST);
            return
                (show_min and show_max) ? output_string(min) + lt_str + type + lt_str + output_string(max) :
                (show_min) ? type + gt_str + output_string(min) :
                (show_max) ? type + lt_str + output_string(max) :
                type;
        }
        virtual bool check(const T &val) const override {
            return strict
                ? (val > min and val < max)
                : (val >= min and val <= max);
        }
        static constexpr T LOWEST = std::numeric_limits<T>::has_infinity ? -std::numeric_limits<T>::infinity() : std::numeric_limits<T>::lowest();
        static constexpr T BIGGEST = std::numeric_limits<T>::has_infinity ? std::numeric_limits<T>::infinity() : std::numeric_limits<T>::max();
        const T min = LOWEST;
        const T max = BIGGEST;
        const bool strict = false;
        const std::string lt_str{strict ? "<" : u8"⩽"}, gt_str{strict ? ">" : u8"⩾"};
};

/** Contains the command line arguments that aren't carried in the Creativity object */
struct cmd_args {
    unsigned int periods, max_threads;
    std::string out;
    bool overwrite;
};
cmd_args parseCmdArgs(int argc, char **argv, Creativity &cr) {
    try {
        CmdLineHack cmd("Eris-based creativity simulation model", ' ', "0.0.1");

        cmd.startReverseHack();

#define ARG(NAME, PARAM, SHORT, LONG, DESC) \
        TCLAP::ValueArg<decltype(cr.parameters.PARAM)> opt_##NAME##_arg(SHORT, LONG, DESC, false, cr.parameters.PARAM, &opt_##NAME##_constr, cmd)
#define OPTION_UNBOUNDED_NAME(NAME, PARAM, SHORT, LONG, DESC) \
        RangeConstraint<decltype(cr.parameters.PARAM)> opt_##NAME##_constr; \
        ARG(NAME,PARAM,SHORT,LONG,DESC)
#define OPTION_BOUND_TYPE(NAME, PARAM, SHORT, LONG, DESC, METHOD, BOUND) \
        auto opt_##NAME##_constr = RangeConstraint<decltype(cr.parameters.PARAM)>::METHOD(BOUND); \
        ARG(NAME,PARAM,SHORT,LONG,DESC)
#define OPTION_LBOUND_NAME(NAME, PARAM, SHORT, LONG, DESC, LBOUND) OPTION_BOUND_TYPE(NAME, PARAM, SHORT, LONG, DESC, GE, LBOUND)
#define OPTION_UBOUND_NAME(NAME, PARAM, SHORT, LONG, DESC, UBOUND) OPTION_BOUND_TYPE(NAME, PARAM, SHORT, LONG, DESC, LE, UBOUND)
#define OPTION_BOUND_NAME(NAME, PARAM, SHORT, LONG, DESC, LBOUND, UBOUND) \
        RangeConstraint<decltype(cr.parameters.PARAM)> opt_##NAME##_constr(LBOUND, UBOUND); \
        ARG(NAME,PARAM,SHORT,LONG,DESC)
#define OPTION_LBOUND(PARAM, SHORT, LONG, DESC, LBOUND) OPTION_LBOUND_NAME(PARAM, PARAM, SHORT, LONG, DESC, LBOUND)
#define OPTION_UBOUND(PARAM, SHORT, LONG, DESC, UBOUND) OPTION_UBOUND_NAME(PARAM, PARAM, SHORT, LONG, DESC, UBOUND)
#define OPTION_LBOUND_STRICT(PARAM, SHORT, LONG, DESC, LBOUND) OPTION_BOUND_TYPE(PARAM, PARAM, SHORT, LONG, DESC, GT, LBOUND)
#define OPTION_UBOUND_STRICT(PARAM, SHORT, LONG, DESC, UBOUND) OPTION_BOUND_TYPE(PARAM, PARAM, SHORT, LONG, DESC, LT, UBOUND)
#define OPTION_BOUND(PARAM, SHORT, LONG, DESC, LBOUND, UBOUND) OPTION_BOUND_NAME(PARAM, PARAM, SHORT, LONG, DESC, LBOUND, UBOUND)
#define OPTION_UNBOUNDED(PARAM, SHORT, LONG, DESC) OPTION_UNBOUNDED_NAME(PARAM, PARAM, SHORT, LONG, DESC)
#define OPTION_LBOUND_INIT(PARAM, SHORT, LONG, DESC, LBOUND) OPTION_LBOUND_NAME(initial_##PARAM, initial.PARAM, SHORT, LONG, DESC, LBOUND)
#define OPTION_UBOUND_INIT(PARAM, SHORT, LONG, DESC, UBOUND) OPTION_UBOUND_NAME(initial_##PARAM, initial.PARAM, SHORT, LONG, DESC, UBOUND)
#define OPTION_BOUND_INIT(PARAM, SHORT, LONG, DESC, LBOUND, UBOUND) OPTION_BOUND_NAME(initial_##PARAM, initial.PARAM, SHORT, LONG, DESC, LBOUND, UBOUND)
#define OPTION_UNBOUNDED_INIT(PARAM, SHORT, LONG, DESC) OPTION_UNBOUNDED_NAME(initial_##PARAM, initial.PARAM, SHORT, LONG, DESC)

// Single-letter options used:
// b B c C d D f i j k K m M n N o O P Q r R s T w W x y z Z

        OPTION_LBOUND(readers, "r", "readers", "Number of reader/author agents in the simulation", 1);
        OPTION_LBOUND(dimensions, "D", "dimensions", "Number of dimensions of the simulation", 1);
        auto opt_density_constr = RangeConstraint<double>::GT(0);
        TCLAP::ValueArg<double> opt_density_arg("d", "density", "Reader density (in readers per unit^(D), where D is the configured # of dimensions)",
                false, Creativity::densityFromBoundary(cr.parameters.readers, cr.parameters.dimensions, cr.parameters.boundary), &opt_density_constr, cmd);
        OPTION_BOUND(piracy_link_proportion, "f", "piracy-link-proportion", "Proportion of potential sharing links between readers that are created", 0, 1);
        OPTION_LBOUND(book_distance_sd, "B", "book-distance-sd", "Standard deviation of book distance from author; distance ~ |N(0, B)|", 0);
        OPTION_LBOUND(book_quality_sd, "Q", "book-quality-sd", "Standard deviation of book perceived quality; perceived quality ~ N(q, Q), where q is the innate quality and Q is this value.", 0);
        OPTION_LBOUND(reader_step_sd, "R", "reader-step-sd", "Standard deviation of the inter-period random-direction reader movement distance ~ |N(0, R)|", 0);
        OPTION_UBOUND_STRICT(reader_creation_shape, "s", "reader-creation-shape", "Shape parameter, β, of the creator effort function: q(ℓ) = (α/β)[(ℓ+1)^β - 1]", 1);
        OPTION_LBOUND(reader_creation_scale_min, "z", "reader-creation-scale-min", "Minimum support of scale parameter α ~ U[a,b] for the creator effort function: q(ℓ) = (α/β)[(ℓ+1)^β - 1]", 0);
        OPTION_LBOUND(reader_creation_scale_max, "Z", "reader-creation-scale-max", "Maximum support of scale parameter α ~ U[a,b] for the creator effort function: q(ℓ) = (α/β)[(ℓ+1)^β - 1]", 0);
        OPTION_LBOUND(cost_fixed, "C", "cost-fixed", "Fixed cost of keeping a book on the market for a period", 0);
        OPTION_LBOUND(cost_unit, "c", "cost-unit", "Unit cost of making a copy of a book", 0);
        OPTION_LBOUND(cost_piracy, "y", "cost-piracy", "Cost of receiving a pirated copy of a book", 0);
        OPTION_LBOUND_STRICT(income, "i", "income", "Per-period external reader income", 0);
        OPTION_BOUND(prior_weight, "w", "prior-weight", "The per-period precision matrix multiplier when using a belief as the next period's prior", 0, 1);
        OPTION_BOUND(prior_weight_piracy, "W", "prior-weight-piracy", "The per-period precision matrix multiplier for the first piracy period", 0, 1);

        OPTION_UNBOUNDED_INIT(belief_threshold, "b", "belief-threshold", "The n-k value at which a readers bases decision on beliefs instead of initial parameters");
        OPTION_BOUND_INIT(prob_write, "x", "initial-prob-write", "The probability of writing in initial periods", 0, 1);
        OPTION_LBOUND_INIT(q_min, "m", "initial-quality-min", "The minimum support of quality q ~ U[a,b] for authored books in initial periods", 0);
        OPTION_LBOUND_INIT(q_max, "M", "initial-quality-max", "The maximum support of quality q ~ U[a,b] for authored books in initial periods", 0);
        OPTION_LBOUND_INIT(p_min, "n", "initial-price-min", "The minimum support of price p ~ U[a,b] for new books in initial periods", 0);
        OPTION_LBOUND_INIT(p_max, "N", "initial-price-max", "The maximum support of price p ~ U[a,b] for new books in initial periods", 0);
        OPTION_BOUND_INIT(prob_keep, "k", "initial-prob-keep", "The probability of keeping a previously-written book on the market for another period", 0, 1);
        OPTION_BOUND_INIT(keep_price, "K", "keep-price", "The  price-above-marginal-cost level (relative to current P-MC) for a book left on the market for another period", 0, 1);

        OPTION_UNBOUNDED(piracy_begins, "P", "piracy-begins", "The period in which piracy becomes available");

        auto periods_constr = RangeConstraint<unsigned int>::GE(0);
        TCLAP::ValueArg<unsigned int> periods_arg("T", "periods", "Number of simulation periods to run", false, 200, &periods_constr, cmd);

        TCLAP::ValueArg<std::string> output_file("o", "output", "Output file (or database URL) for simulation results.  Example database URL: 'postgresql://user:secret@localhost:5432/dbname?sslmode=require&otheroption=123'", false,
                "creativity-" + std::to_string(Random::seed()) + ".crstate", // default
                "filename-or-database", cmd);

        TCLAP::SwitchArg overwrite_output_file("O", "overwrite", "Allows output file given to -o to be overwritten.  No effect if a database URL is given to -o.", cmd, false);

        RangeConstraint<unsigned int> max_threads_constr(0, std::thread::hardware_concurrency());
        TCLAP::ValueArg<unsigned int> max_threads_arg("j", "threads", "Maximum number of threads to use for the simulation", false, 0, &max_threads_constr, cmd);

        cmd.reverseHack();

        cmd.parse(argc, argv);

#define COPY_PARAM_SETTING(PARAM) cr.set().PARAM = opt_##PARAM##_arg.getValue()
#define COPY_INIT_PARAM_SETTING(PARAM) cr.set().initial.PARAM = opt_initial_##PARAM##_arg.getValue()
        COPY_PARAM_SETTING(readers);
        COPY_PARAM_SETTING(dimensions);
        cr.set().boundary = Creativity::boundaryFromDensity(cr.parameters.readers, cr.parameters.dimensions, opt_density_arg.getValue());
        COPY_PARAM_SETTING(book_distance_sd);
        COPY_PARAM_SETTING(book_quality_sd);
        COPY_PARAM_SETTING(reader_step_sd);
        COPY_PARAM_SETTING(reader_creation_shape);
        COPY_PARAM_SETTING(reader_creation_scale_min);
        COPY_PARAM_SETTING(reader_creation_scale_max);
        COPY_PARAM_SETTING(cost_fixed);
        COPY_PARAM_SETTING(cost_unit);
        COPY_PARAM_SETTING(cost_piracy);
        COPY_PARAM_SETTING(income);
        COPY_PARAM_SETTING(piracy_begins);
        COPY_PARAM_SETTING(piracy_link_proportion);
        COPY_PARAM_SETTING(prior_weight);
        COPY_PARAM_SETTING(prior_weight_piracy);
        COPY_INIT_PARAM_SETTING(prob_write);
        COPY_INIT_PARAM_SETTING(q_min);
        COPY_INIT_PARAM_SETTING(q_max);
        COPY_INIT_PARAM_SETTING(p_min);
        COPY_INIT_PARAM_SETTING(p_max);
        COPY_INIT_PARAM_SETTING(prob_keep);
        COPY_INIT_PARAM_SETTING(keep_price);
        COPY_INIT_PARAM_SETTING(belief_threshold);

        cmd_args ret;
        ret.periods = periods_arg.getValue();
        ret.max_threads = max_threads_arg.getValue();
        ret.out = output_file.getValue();
        ret.overwrite = overwrite_output_file.getValue();

        return ret;
    }
    catch (TCLAP::ArgException &e) {
        std::cerr << "Error: " << e.error() << " for arg " << e.argId() << "\n";
        throw;
    }
}


int main(int argc, char *argv[]) {
    std::cerr << std::setprecision(16);
    std::cout << std::setprecision(16);
    Eigen::initParallel();
    auto creativity = Creativity::create();

    auto args = parseCmdArgs(argc, argv, *creativity);

    std::string results_out_display = args.out;

    if (args.out.substr(0, 13) == "postgresql://" or args.out.substr(0, 11) == "postgres://") {
        try {
            creativity->pgsql(args.out, false /*read-only*/, true /*write-only*/);

#ifndef CREATIVITY_SKIP_PGSQL
            PsqlStorage &pgsql = dynamic_cast<PsqlStorage&>(creativity->storage().first->backend());
            auto conn_locked = pgsql.connection();
            auto &conn = conn_locked.first;
            // Rebuild an output URL for the postgresql database
            std::ostringstream psqlurl("postgresql://");
            const char *username = conn.username();
            if (username) psqlurl << username << "@";
            // Omit password
            const char *hostname = conn.hostname();
            if (hostname) psqlurl << hostname;
            const char *port = conn.port();
            if (port) psqlurl << ":" << port;
            psqlurl << "/" << conn.dbname() << "?creativity=" << pgsql.id;
            results_out_display = psqlurl.str();
#endif
        }
        catch (const std::exception &e) {
            std::cerr << "Connection to PostgreSQL failed: " << e.what() << "\n";
            exit(1);
        }
    }
    else {
        // Filename output

        try {
            // If the user didn't specify --overwrite, make sure the file doesn't exist
            if (not args.overwrite) {
                struct stat buffer;
                if (stat(args.out.c_str(), &buffer) == 0) { // The file exists
                    std::cerr << "Error: `" << args.out << "' already exists; specify a different file or add `--overwrite' option to overwrite\n";
                    exit(1);
                }
            }
            creativity->fileWrite(args.out);
        }
        catch (std::exception &e) {
            std::cerr << "Unable to write to file: " << e.what() << "\n";
            exit(1);
        }
    }
    std::cout << "Writing simulation results to " << results_out_display << "\n";

    creativity->setup();
    auto sim = creativity->sim;

    // Copy the initial state into the storage object
    creativity->storage().first->emplace_back(sim);

    sim->maxThreads(args.max_threads);

    std::chrono::time_point<std::chrono::high_resolution_clock> now,
        last = std::chrono::high_resolution_clock::now();

    while (sim->t() < args.periods) {
        sim->run();

        creativity->storage().first->emplace_back(sim);

        now = std::chrono::high_resolution_clock::now();
        double speed = 1.0 / std::chrono::duration<double>(now - last).count();
        std::swap(last, now);
        std::cout << "\rRunning simulation [t=" << sim->t() << "; R=" << sim->countAgents<Reader>() << "; B=" << sim->countGoods<Book>() << ", newB=" <<
            sim->countGoods<Book>([](const Book &b) -> bool { return b.age() == 0; }) << "] " << speed << " Hz                    " << std::flush;
    }

    creativity->storage().first->flush();

    std::cout << "\n\nSimulation complete.  Results saved to " << results_out_display << "\n\n";
}
