#pragma once
#include "creativity/cmdargs/Validation.hpp"
#include "creativity/cmdargs/strings.hpp"
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <string>
#include <limits>
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
        void parse(int argc, char const* const* argv);

        /** Creates an option value object without any special validation wrapper class.  This
         * function participates only when `T` is not an unsigned type.
         *
         * \param store the default value and the location to store a specified value.
         */
        template <typename T>
        static typename std::enable_if<not std::is_unsigned<T>::value and not std::is_same<T, bool>::value,
                                       boost::program_options::typed_value<T>*
                                      >::type
        value(T& store) {
            return boost::program_options::value<T>(&store)->default_value(store)->value_name(type_string<T>());
        }

        /** Creates an option value object around an unsigned primitive type, with automatic value
         * storage and default value.  This function participates only when `T` is an unsigned type,
         * and, functionally, will ensure that a negative value such as -3 is not accepted, but will not
         * otherwise restrict the value.
         */
        template <typename T>
        static typename std::enable_if<std::is_unsigned<T>::value and not std::is_same<T, bool>::value,
                        boost::program_options::typed_value<Validation<T>>*>::type
        value(T &storage) {
            return value<Validation<T>>(storage);
        }

        /** Creates an option value for a boolean value, that is, for an switch without an argument,
         * with default value as given in `store`.
         */
        static boost::program_options::typed_value<bool>* value(bool& store) {
            return boost::program_options::bool_switch(&store)->default_value(store);
        }

        /** Creates an option value object with explicit validation wrapper class V.  `store` is
         * used for both the default and the location to store a command-line provided value.
         *
         * \tparam V a class that will throw during construction if the given value isn't valid; the
         * class must expose a value type in `value_type`, which should match the given `store`
         * value.
         * \param store the location of the default value and the location to store a given value
         */
        template <typename V>
        static boost::program_options::typed_value<V>*
        value(typename V::value_type &store) {
            return boost::program_options::value<V>()->default_value(store)->value_name(V::validationString())
                ->notifier([&store](const V &v) { store = v; /* NB: implicit conversion */ });
        }

        /** Creates an option value object around a vector of options with validation wrapper class
         * V applied to each element of the vector.  `store` is used both for the default values and
         * the location to store command-line provided values.
         *
         * \tparam V a class that will throw during construction if the given value isn't valid; the
         * class must expose a value type in `value_type`, which must match the `value_type` of the
         * given vector.
         * \param store the vector in which to store given values.
         */
        template <typename V, typename A>
        static boost::program_options::typed_value<std::vector<V>>*
        value(std::vector<typename V::value_type, A> &store) {
            return boost::program_options::value<std::vector<V>>()->value_name(
                    V::validationString() + " [" + V::validationString() + " ...]")
                ->notifier([&store](const std::vector<V> &v) {
                        store.clear();
                        store.reserve(v.size());
                        for (const auto &val : v) store.push_back((const typename V::value_type) val);
                    });
        }

        /** Takes a std::vector of values for options that store multiple values. This function
         * participates only when the vector stores non-unsigned types. */
        template <typename T, typename A>
        typename std::enable_if<not std::is_unsigned<T>::value, boost::program_options::typed_value<std::vector<T, A>>*>::type
        value(std::vector<T, A> &store) {
            return boost::program_options::value<std::vector<T, A>>(&store)->value_name(
                    type_string<T>() + " [" + type_string<T>() + " ...]");
        }

        /** Takes a std::vector of values for options that store multiple values. This function
         * participates only when the vector stores unsigned types. */
        template <typename T, typename A>
        typename std::enable_if<std::is_unsigned<T>::value, boost::program_options::typed_value<std::vector<Validation<T>>>*>::type
        value(std::vector<T, A> &store) {
            return value<Validation<T>>(store);
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

        /** Returns a usage string such as "Usage: program [ARGS]".  Called by help().  Subclasses should
         * override to change the string as needed.
         */
        virtual std::string usage() const;

        /// Returns a argument help message.
        virtual std::string help() const;

    protected:
        /** Adds options.  This method is called automatically by parse() before parsing arguments
         * if nothing has been set in the options_ object; the default implementation adds --help
         * and --version options.  Subclasses should override and enhance to also populate
         * `options_` appropriately.
         */
        virtual void addOptions();

        /** The program name, populated by parse(). */
        std::string prog_name_;

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





}}
