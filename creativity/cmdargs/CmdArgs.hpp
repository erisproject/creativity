#pragma once
#include "creativity/Creativity.hpp"
#include "creativity/cmdargs/Validation.hpp"
#include "creativity/cmdargs/strings.hpp"
#include <boost/program_options/errors.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <eris/Random.hpp>
#include <string>
#include <limits>
#include <regex>
#include <type_traits>
#include <vector>

namespace boost { namespace program_options { class variables_map; } }

namespace creativity {
/// Namespace for command-line argument handling classes.
namespace cmdargs {

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

    protected:
        /// Default constructor is protected; use a suitable subclass.
        CmdArgs() = default;

    public:
        /** Parses options and verifies them.  Throws an std::exception-derived exception for
         * various errors (such as invalid arguments or invalid argument types).
         *
         * A suitable method such as addCliOptions() or addGuiOptions() must have been called first:
         * otherwise there will be no recognized options to parse.
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
            return boost::program_options::value<T>(&store)->default_value(store)->value_name(type_string<T>());
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

        /// Returns a version string.
        virtual std::string version() const;

        /// Returns a argument help message.
        virtual std::string help() const;

    protected:
        /** Adds options.  This method is called automatically by parse() before parsing arguments
         * if nothing has been set in the options_ object; the default implementation adds --help
         * and --version options.  Subclasses should override to also populate `options_`
         * appropriately.
         */
        virtual void addOptions();

        /// The options descriptions variable for all options.
        boost::program_options::options_description options_;

        /// Like options_, but for hidden options that aren't to be displayed in --help output.
        boost::program_options::options_description invisible_;

        /// Positional options object
        boost::program_options::positional_options_description positional_;

        /** Does nothing; subclasses should override if needed to deal with argument values.  This
         * is called after checking for and handling --help or --version flags.
         *
         * \param vars the parsed variables
         */
        virtual void postParse(boost::program_options::variables_map &vars);
};


/** Overload of validate for boost to convert from string to a validated data type.  Mostly this
 * just checks that the validation object can be constructed (which will fail if the validation
 * fails), but this also makes sure signed types aren't provided with a leading minus sign.
 */
template <class V, typename = typename std::enable_if<std::is_base_of<ValidationTag, V>::value>::type>
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

}}
