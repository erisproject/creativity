#include "creativity/Creativity.hpp"
#include <eris/Simulation.hpp>
#include <eris/Random.hpp>
#include <functional>
#include <iostream>
#include <iomanip>
#include <regex>
#include <Eigen/Core>
#include <tclap/CmdLine.h>
#include <sys/stat.h>

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
        RangeConstraint() = delete;
        RangeConstraint(T min) : min_(min), no_max_(true) {}
        RangeConstraint(T min, T max) : min_(min), max_(max) {}
        virtual std::string description() const override {
            if (max_ == std::numeric_limits<double>::infinity())
                return "Value must be at least " + output_string(min_);
            else
                return "Value must be between " + output_string(min_) + " and " + output_string(max_);
        }
        virtual std::string shortID() const override {
            if (no_max_)
                return "V ≥ " + output_string(min_);
            else
                return "V ∈ [" + output_string(min_) + "," + output_string(max_) + "]";
        }
        virtual bool check(const T &val) const override {
            return val >= min_ and val <= max_;
        }
    private:
        T min_;
        T max_ = std::numeric_limits<T>::has_infinity ? std::numeric_limits<T>::infinity() : std::numeric_limits<T>::max();
        bool no_max_ = false;
};

/// Constraint for an array of a given length.
class ArrayConstraint : public TCLAP::Constraint<std::string> {
    public:
        ArrayConstraint(unsigned int size) : size_(size) {
            if (size < 2) throw std::logic_error("ArrayConstraint(s) requires s >= 2");
        }
        virtual std::string description() const override {
            return "Value must be a comma-separated list of " + std::to_string(size_) + " numeric values";
        }
        virtual std::string shortID() const override {
            return "NUMERIC[" + std::to_string(size_) + "]";
        }
        virtual bool check(const std::string &val) const override {
            std::regex re("^(?:(?:^|,)" + number_re_ + "){" + std::to_string(size_) + "}$");
            return std::regex_match(val, re);
        }
    private:
        const std::string number_re_{"(?:\\d+(?:\\.\\d+)?|\\.\\d+)(?:[eE][-+]?\\d+)?"};
        const unsigned int size_;
};

/// Returns vector elements concatenated together with a ,
std::string vectorString(const Ref<const VectorXd> &v) {
    std::ostringstream str;
    for (int i = 0; i < v.cols(); i++) {
        if (i > 0) str << ",";
        str << v[i];
    }
    return str.str();
}

/// Returns lower triangle matrix elements in row-major order
std::string lowerTriangleString(const Ref<const MatrixXd> &M) {
    if (M.cols() != M.rows()) throw std::logic_error("lowerTriangleString requires a square matrix");
    std::ostringstream str;
    for (int r = 0; r < M.rows(); r++) {
        for (int c = 0; c <= r; c++) {
            if (r > 0 or c > 0) str << ",";
            str << M(r,c);
        }
    }
    return str.str();
}

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

#define OPTION_LBOUND(PARAM, SHORT, LONG, DESC, LBOUND) \
        RangeConstraint<decltype(cr.parameters.PARAM)> opt_##PARAM##_constr(LBOUND); \
        TCLAP::ValueArg<decltype(cr.parameters.PARAM)> opt_##PARAM##_arg(SHORT, LONG, DESC, false, cr.parameters.PARAM, &opt_##PARAM##_constr, cmd);

        OPTION_LBOUND(dimensions, "D", "dimensions", "Number of dimensions of the simulation", 1);
        OPTION_LBOUND(readers, "r", "readers", "Number of reader/author agents in the simulation", 2);
        OPTION_LBOUND(density, "d", "density", "Reader density (in readers per unit^(D), where D is the configured # of dimensions)",
                std::numeric_limits<double>::min());
        OPTION_LBOUND(book_distance_sd, "B", "book-distance-sd", "Standard deviation of book distance from author; distance ~ |N(0, B)|", 0);
        OPTION_LBOUND(book_quality_sd, "Q", "book-quality-sd", "Standard deviation of book perceived quality; perceived quality ~ N(q, Q), where q is the innate quality and Q is this value.", 0);
        OPTION_LBOUND(cost_fixed, "C", "cost-fixed", "Fixed cost of keeping a book on the market for a period", 0);
        OPTION_LBOUND(cost_unit, "c", "cost-unit", "Unit cost of making a copy of a book", 0);
        OPTION_LBOUND(income, "i", "income", "Per-period external reader income", std::numeric_limits<double>::min());

#define OPTION_VECTOR(PARAM, SHORT, LONG, DESC) \
        ArrayConstraint opt_##PARAM##_constr(cr.parameters.PARAM.rows()); \
        TCLAP::ValueArg<std::string> opt_##PARAM##_arg(SHORT, LONG, DESC, false, vectorString(cr.parameters.PARAM), &opt_##PARAM##_constr, cmd);
#define OPTION_LOWERTRI(PARAM, SHORT, LONG, DESC) \
        ArrayConstraint opt_##PARAM##_constr(cr.parameters.PARAM.rows() * (cr.parameters.PARAM.rows() + 1) / 2); \
        TCLAP::ValueArg<std::string> opt_##PARAM##_arg(SHORT, LONG, DESC, false, lowerTriangleString(cr.parameters.PARAM), &opt_##PARAM##_constr, cmd);

        RangeConstraint<unsigned int> periods_constr(1);
        TCLAP::ValueArg<unsigned int> periods_arg("p", "periods", "Number of simulation periods to run", false, 100, &periods_constr, cmd);

        TCLAP::ValueArg<std::string> output_file("o", "output", "Output file for simulation results", false,
                "creativity-" + std::to_string(Random::seed()) + ".crstate", // default
                "filename", cmd);

        TCLAP::SwitchArg overwrite_output_file("O", "overwrite", "Allows output file given to -o to be overwritten.", cmd, false);

        RangeConstraint<unsigned int> max_threads_constr(0, std::thread::hardware_concurrency());
        TCLAP::ValueArg<unsigned int> max_threads_arg("T", "threads", "Maximum number of threads to use for the simulation", false, 0, &max_threads_constr, cmd);

        cmd.reverseHack();

        cmd.parse(argc, argv);

#define COPY_PARAM_SETTING(PARAM) cr.parameters.PARAM = opt_##PARAM##_arg.getValue()
        COPY_PARAM_SETTING(dimensions);
        COPY_PARAM_SETTING(readers);
        COPY_PARAM_SETTING(density);
        COPY_PARAM_SETTING(book_distance_sd);
        COPY_PARAM_SETTING(book_quality_sd);
        COPY_PARAM_SETTING(cost_fixed);
        COPY_PARAM_SETTING(cost_unit);
        COPY_PARAM_SETTING(income);

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
    // If the user didn't specify --overwrite, make sure the file doesn't exist
    if (not args.overwrite) {
        struct stat buffer;
        if (stat(args.out.c_str(), &buffer) == 0) { // The file exists
            std::cerr << "Error: `" << args.out << "' already exists; specify a different file or add `--overwrite' option to overwrite\n";
        }
    }
    creativity->fileWrite(args.out);

    std::cout << "Writing simulation results to `" << args.out << "'\n";

    creativity->setup();
    auto sim = creativity->sim;

    // Copy the initial state into the storage object
    creativity->storage().first->emplace_back(sim);

    ERIS_DBG("here we go...!");

    while (sim->t() < args.periods) {
        ERIS_DBGVAR(sim->t());
        ERIS_DBGVAR(args.periods);
        ERIS_DBG("running");
        sim->run();
        ERIS_DBG("done running");

        creativity->storage().first->emplace_back(sim);
        ERIS_DBGVAR(sim->t());
    }
    ERIS_DBGVAR(creativity->storage().first->size());

    creativity->storage().first->flush();

    std::cout << "\n\nSimulation complete.  Results saved to " << args.out << "\n\n";
}
