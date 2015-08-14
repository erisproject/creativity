#pragma once
#include <boost/program_options.hpp>
#include <list>
#include <limits>
#include <regex>
#include <eris/Random.hpp>
#include "creativity/Creativity.hpp"

namespace creativity {

/// constexpr that returns positive infinity, if T has such a value, or the maximum value T supports if it does not.
template <typename T> constexpr T inf_or_max() { return std::numeric_limits<T>::has_infinity ? std::numeric_limits<T>::infinity() : std::numeric_limits<T>::max(); }
/// constexpr that returns negative infinity, if T has such a value, or the lowest value T supports if it does not.
template <typename T> constexpr T neginf_or_lowest() { return std::numeric_limits<T>::has_infinity ? -std::numeric_limits<T>::infinity() : std::numeric_limits<T>::lowest(); }

/** This class handles command line argument parsing.  Most options are supported both by the cli
 * interface (to specify simulation settings) and by the gui (to specify simulation defaults).
 *
 * \todo Argument output isn't quite right because boost::p_a doesn't know about the utf8 strings
 * being used, so alignment fails.  Figure out a way to fix that, perhaps by using std::wstring.
 *
 * \todo Add settings *plus* colour setting arguments to a file (~/.creativity-settings.rc).
 */
class CmdArgs {

    public:
        /// No default constructor.
        CmdArgs() = delete;

        /// Constructs a CmdArgs object that will use and adjust the given CreativitySettings.
        CmdArgs(CreativitySettings &s) : s_{s} {}

        /** Returns a value converted to a string; in most cases this just passes the value to
         * std::to_string for conversion.
         * 
         * A specialization of this for doubles converts the double to a string via
         * std::to_string(), then trims off trailing 0's and, possibly, a trailing `.`.
         */
        template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
        static std::string output_string(T v);

        /** Parses options and verifies them.  Throws an std::exception-derived exception for
         * various errors (such as invalid arguments or invalid argument types).
         *
         * Either addCliOptions() or addGuiOptions() must have been called first: otherwise there
         * will be no recognized options to parse.
         *
         * Options that can be stored directly in the CreativitySettings object passed to the
         * constructor are assigned directly; other variable are accessible in the `variables`
         * member.
         *
         * The --help and --version arguments are handled internally: they print the requested info,
         * then exit the program.
         *
         * \param argc the argc as received by main
         * \param argv the argv as received by main
         */
        void parse(int argc, char *argv[]);

        /** These are the variables that can be set by a call to parse().  Whether or not they are
         * actually set depends on the command-line options made available: for example, `.start` is
         * only set if addGuiOptions() was called, while `.tmpdir` is only for addCliOptions().
         *
         * Note that simulation-specific settings are set directly into the CreativitySettings
         * object given during construction.
         */
        struct {
            /** The number of periods to run the simulation.  Before calling addCliOptions() or
             * addGuiOptions(), this is the default value; after calling parse this will be updated
             * to whatever the user specified.
             */
            unsigned int periods = 250;
            /// Whether to start running the simulation in the GUI right away
            bool start = false;
            /** Whether to initialize (but not start) the simulation in the GUI right away.  Has no
             * effect if start is true.
             */
            bool initialize = false;
            /// The output file for simulation results; has a default for the CLI, not for the GUI
            std::string output;
            /// The temporary directory for results (supported by CLI only)
            std::string tmpdir;
            /// Whether `output' can be overwritten (supported by CLI only)
            bool overwrite = false;
            /// The input file (for the GUI)
            std::string input;
            /** The number of threads to use.  The default is 0 for CLI, number of CPU threads for
             * the GUI.
             */
            unsigned int threads = 0;
            /** The seed.  The default is whatever eris::Random::seed() returns, which is random
             * (unless overridden with ERIS_RNG_SEED).  This value can be ignored: it is handled by
             * parse().
             */
            typename eris::Random::rng_t::result_type seed = eris::Random::seed();
        } parameters;

        /// Adds CLI command-line options into the option descriptions
        void addCliOptions();

        /// Adds GUI command-line options into the option descriptions
        void addGuiOptions();

        /** Validation tag; any Validation class must (ultimately) inherit from this class.  This
         * class does nothing.
         */
        class ValidationTag {};

        /** Validation wrapper base class.  This base class does no actual value validation (except
         * one, see below).
         *
         * This should be inherited from virtually so that subclasses can inherit from multiple
         * Validation subclasses to enforce multiple validations at once.
         *
         * Note that there is one validation that occurs when using this: if `T` is an unsigned
         * type, validation of a Validation (or derived class) makes sure the given value is not
         * negative; without this, boost will cast negative input values to unsigned types, thus
         * ending up with -1 becoming (for an unsigned int) 4294967295.
         */
        template <typename T>
        class Validation : public ValidationTag {
            public:
                /// Constructs with an initial value
                Validation(T v) : val_(v) {}
                /// Implicit conversion to the stored value
                operator const T& () const { return val_; }
                /// The type T that this object validates
                using value_type = T;
                /// Returns a string representation of this validation object
                static std::string validationString() {
                    if (std::is_unsigned<T>::value) { return typeString<T>() + u8"⩾0"; }
                    return typeString<T>();
                }
                /// Virtual destructor
                virtual ~Validation() = default;
            protected:
                /// The stored value
                T val_;
        };

        /** Validation wrapper for options that have a minimum value.
         *
         * For technical reasons, the minimum must be given as a fraction of longs (non-type
         * template parameters cannot be non-integral types).
         *
         * \param T any numeric type.
         * \param min the minimum accepted value; if the value is fractional, this is the
         * numerator.
         * \param denom the denominator of the minimum value; defaults to 1.
         */
        template <typename T, long min, long denom = 1>
        class Min : public virtual Validation<T> {
            public:
                /// Constructor.  Throws if `v < min`.
                Min(T v) : Validation<T>(v) {
                    if (*this < min / (double) denom)
                        throw boost::program_options::validation_error(boost::program_options::validation_error::invalid_option_value);
                }

                /// Returns string representation of this validation
                static std::string validationString() { return typeString<T>() + u8"⩾" + (denom == 1 ? output_string(min) : output_string(min / (double) denom)); }
        };

        /** Validation wrapper for options that have a maximum value.
         *
         * For technical reasons, the minimum must be given as a fraction of longs (non-type
         * template parameters cannot be non-integral types).
         *
         * \param T any numeric type.
         * \param max the maximum accepted value; if fractional, this is the numerator.
         * \param denom the denominator of the maximum value; defaults to 1.
         */
        template<typename T, long max, long denom = 1>
        class Max : public virtual Validation<T> {
            public:
                /// Constructor.  Throws if `v > max`.
                Max(T v) : Validation<T>(v) {
                    if (*this > max / (double) denom)
                        throw boost::program_options::validation_error(boost::program_options::validation_error::invalid_option_value);
                }

                /// Returns string representation of this validation
                static std::string validationString() { return typeString<T>() + u8"⩽" + (denom == 1 ? output_string(max) : output_string(max / (double) denom)); }
        };

        /** Validation wrapper for options that have both a minimum and maximum value.
         *
         * \param T any numeric type
         * \param min the minimum accepted value; if fractional, this is the numerator.
         * \param max the maximum accepted value; if fractional, this is the numerator.
         * \param denom the denominator of both the min and max values; defaults to 1.
         */
        template<typename T, long min, long max, long denom = 1>
        class Range : public Min<T, min>, public Max<T, max> {
            public:
                /// Constructor.  Throws if `v < min` or `v > max`.
                Range(T v) : Validation<T>(v), Min<T, min, denom>(v), Max<T, max, denom>(v) {}

                /// Returns string representation of this validation
                static std::string validationString() { return (denom == 1 ? output_string(min) : output_string(min / (double) denom)) + u8"⩽" + typeString<T>() + u8"⩽" + (denom == 1 ? output_string(max) : output_string(max / (double) denom)); }
        };

        /** Validation wrapper for options that have a strict inequality boundary that the value
         * must be above, such as \f$v > 3\f$.
         *
         * \param T any numeric type
         * \param lower the lower bound that the value must be above (or the numerator of the value)
         * \param the denominator under `lower`.  Defaults to 1.
         */
        template <typename T, long lower, long denom = 1>
        class Above : public virtual Validation<T> {
            public:
                /// Constructor.  Throws if `v <= lower`.
                Above(T v) : Validation<T>(v) {
                    if (*this <= lower / (double) denom)
                        throw boost::program_options::validation_error(boost::program_options::validation_error::invalid_option_value);
                }

                /// Returns string representation of this validation
                static std::string validationString() { return typeString<T>() + u8">" + (denom == 1 ? output_string(lower) : output_string(lower / (double) denom)); }
        };

        /** Validation wrapper for options that have a strict inequality boundary that the value
         * must be below, such as \f$v < 3\f$.
         *
         * \param T any numeric type
         * \param upper the upper bound that the value must be above (or the numerator of the value)
         * \param the denominator under `lower`.  Defaults to 1.
         */
        template <typename T, long upper, long denom = 1>
        class Below : public virtual Validation<T> {
            public:
                /// Constructor.  Throws if `v >= upper`.
                Below(T v) : Validation<T>(v) {
                    if (*this >= upper / (double) denom)
                        throw boost::program_options::validation_error(boost::program_options::validation_error::invalid_option_value);
                }

                /// Returns string representation of this validation
                static std::string validationString() { return typeString<T>() + u8"<" + (denom == 1 ? output_string(upper) : output_string(upper / (double) denom)); }
        };

        /** Creates an option value object with validation wrapper class V.  `store` is used for
         * both the default and the location to store a given value.
         *
         * \tparam V a class that will throw during construction if the given value isn't valid; the
         * class must expose a value type in `value_type`, which should match the given `store`
         * value.
         * \param store the location of the default value and the location to store a given value
         */
        template <typename V>
        static boost::program_options::typed_value<V>* value(typename V::value_type& store) {
            return boost::program_options::value<V>()->default_value(store)->value_name(V::validationString())
                ->notifier([&store](const V &v) { store = v; /* NB: implicit conversion */ });
        }

        /** Creates an option value object without any special validation wrapper class.  This
         * function participates only when `T` is not an unsigned type.
         *
         * \param store the default value and the location to store a specified value.
         */
        template <typename T>
        static typename std::enable_if<not std::is_base_of<ValidationTag, T>::value and not std::is_unsigned<T>::value,
                        boost::program_options::typed_value<T>*>::type
        value(T& store) {
            return boost::program_options::value<T>(&store)->default_value(store)->value_name(typeString<T>());
        }

        /** Creates an option value object around an unsigned primitive type, with automatic value
         * storage and default value.  This function participates only when `T` is an unsigned type,
         * and, functionally, will ensure that a value such as -3 is not accepted, but will not
         * otherwise restrict the value.
         */
        template <typename T>
        static typename std::enable_if<not std::is_base_of<ValidationTag, T>::value and std::is_unsigned<T>::value,
                        boost::program_options::typed_value<Validation<T>>*>::type
        value(T &storage) {
            return value<Validation<T>>(storage);
        }

        /// Shortcut for `value<Min<T, n, d>>(val)` with `T` last (so that it can be inferred from `val`)
        template <long minimum, long denom = 1, typename T>
        static inline boost::program_options::typed_value<Min<T, minimum, denom>>* min(T &store) { return value<Min<T, minimum, denom>>(store); }

        /// Shortcut for `value<Max<T, n, d>>(val)` with `T` last (so that it can be inferred from `val`)
        template <long maximum, long denom = 1, typename T>
        static boost::program_options::typed_value<Max<T, maximum, denom>>* max(T &store) { return value<Max<T, maximum, denom>>(store); }

        /// Shortcut for `value<Range<T, a, b, d>>(val)` with `T` last (so that it can be inferred from `val`)
        template <long min, long max, long denom = 1, typename T>
        static boost::program_options::typed_value<Range<T, min, max, denom>>* range(T &store) { return value<Range<T, min, max, denom>>(store); }

        /// Shortcut for `value<Above<T, a, d>>(val)` with `T` last (so that it can be inferred from `val`)
        template <long lower, long denom = 1, typename T>
        static boost::program_options::typed_value<Above<T, lower, denom>>* above(T &store) { return value<Above<T, lower, denom>>(store); }

        /// Shortcut for `value<Below<T, b, d>>(val)` with `T` last (so that it can be inferred from `val`)
        template <long upper, long denom = 1, typename T>
        static boost::program_options::typed_value<Below<T, upper, denom>>* below(T &store) { return value<Below<T, upper, denom>>(store); }

    protected:

        /** Adds common options into the options descriptions.  Called by
         * addCliOptions()/addGuiOptions().
         */
        void addCommonOptions();

        /// The options descriptions variable for visible options
        boost::program_options::options_description desc_;

        /** The options (which will be added to the end of desc_) containing simulator settings,
         * some of which differ between the CLI and the GUI.
         */
        boost::program_options::options_description sim_desc_{"Simulator settings"};

        /// The options descriptions variable for invisible options
        boost::program_options::options_description invisible_;

        /// Positional options
        boost::program_options::positional_options_description pos_;

        /** Returns an argument name for T.  For numeric types, this is one of 𝑹 (for floating point
         * values), 𝑵 (for unsigned integer types), or 𝒁 (for signed integer types).  For other
         * types, this is simply "arg".
         */
        template <typename T> static std::string typeString();


    private:
        CreativitySettings &s_;

        // Density pseudo-parameter: the stored value is actually boundary
        double density_ = Creativity::densityFromBoundary(s_.readers, s_.dimensions, s_.boundary);
};


template <typename T, typename>
std::string CmdArgs::output_string(T v) {
    return std::to_string(v);
}


/** Specialization of output_string for double types that trims trailing 0's and a trailing decimal
 * point from the string representation before returning it.
 */
template <> std::string CmdArgs::output_string(double v); // Specialization in .cpp

template <typename T> std::string CmdArgs::typeString() {
    return std::is_floating_point<T>::value ? u8"ℝ" :
        std::is_integral<T>::value ? std::is_unsigned<T>::value ? u8"ℕ" : u8"ℤ" :
        u8"arg";
}

/** Overload of validate for boost to convert from string to a validated data type.  Mostly this
 * just checks that the validation object can be constructed (which will fail if the validation
 * fails), but this also makes sure signed types aren't provided with a leading minus sign.
 */
template <class V, typename = typename std::enable_if<std::is_base_of<CmdArgs::ValidationTag, V>::value>::type>
void validate(boost::any &v, const std::vector<std::string> &values, V*, int) {
    using namespace boost::program_options;

    // Check that this value hasn't already been assigned:
    validators::check_first_occurrence(v);
    // Check that only a single value was given, and get it:
    std::string s(validators::get_single_string(values));

    if (std::is_unsigned<typename V::value_type>::value and std::regex_search(s, std::regex("^\\s*-")))
        throw invalid_option_value(s);

    try {
        // First convert to the appropriate type (this will throw if that can't be done):
        auto val = boost::lexical_cast<typename V::value_type>(s);

        // The constructor here will throw if validation fails:
        v = boost::any(V(val));
    }
    catch (...) {
        throw invalid_option_value(s);
    }
}

}

