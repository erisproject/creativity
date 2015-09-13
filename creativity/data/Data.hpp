#pragma once
#include <string>
#include <functional>
#include <vector>
#include <algorithm>
#include <type_traits>
#include <eris/types.hpp>

/** \file
 * Various routines for generating data from stored creativity results.
 */

namespace creativity { struct CreativitySettings; }
namespace creativity { namespace state { class Storage; } }

namespace creativity { namespace data {

/** Struct for holding a simulator parameter. */
struct initial_datum {
    /// The name of this parameter (to be included in the output header).
    const std::string name;
    /// Returns an integer value for this parameter, if callable.
    const std::function<long(const CreativitySettings&)> calc_int;
    /// Returns a double value for this parameter, if callable.
    const std::function<double(const CreativitySettings&)> calc_double;

    /// Default constructor disabled
    initial_datum() = delete;
    /// Constructs a simulation parameter datum that returns an integral value.
    template <typename IntType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
    initial_datum(std::string &&name, std::function<IntType(const CreativitySettings&)> &&calc_int)
        : name(std::move(name)), calc_int(std::move(calc_int)) {}
    /// Constructs a simulation parameter datum that returns a double value.
    initial_datum(std::string &&name, std::function<double(const CreativitySettings&)> &&calc_double)
        : name(std::move(name)), calc_double(std::move(calc_double)) {}
};

/** Struct for holding a simulation data value calculator. */
struct datum {
    /// The name of this parameter (will be prefixed with pre_, new_, or post_)
    const std::string name;
    /** Callable object called to calculate the value of the parameter.
     *
     * The object receives the simulation storage object and time periods `from` and `to`: it should
     * calculate the statistic for the period from `from` to `to`, inclusive.
     */
    const std::function<double(const state::Storage&, eris::eris_time_t, eris::eris_time_t)> calculate;
    /// If true, this statistic is included in the pre-piracy data
    struct {
        bool pre; ///< True if this field applies to pre-piracy periods
        bool piracy; ///< True if this field applies to piracy periods
        bool public_sharing; ///< True if this field applies to public sharing periods
    } applies_to;

    /// Default constructor deleted
    datum() = delete;
    /** Constructor.
     *
     * \param name the parameter name, which will be prefixed with "pre_", "new_", or "post_" when
     * included in the output data.  "pre_" is for the pre-piracy period; "new_" is for the
     * periods just after piracy is introduced; and "post_" is for the final simulation periods.
     * \param calc the callable object to store in `calculate`.
     * \param pre specifies whether this data field applies to pre-piracy periods (default true)
     * \param piracy specifies whether this data field applies to piracy periods (default true)
     * \param public_sharing specifies whether this data field applies to public sharing periods (default true)
     */
    datum(std::string &&name, decltype(calculate) calc, bool pre = true, bool piracy = true, bool public_sharing = true)
        : name(std::move(name)), calculate(std::move(calc)), applies_to{pre, piracy, public_sharing} {}
};

/** Calculates the average net utility over the given period.  Net utility is a reader's utility
 * minus the reader's income (since without any simulation activity, the reader's utility is
 * quasilinear, thus the reader receives exactly his income as utility).
 */
double net_u(const state::Storage &cs, eris::eris_time_t from, eris::eris_time_t to);

/** Calculates the average market life of books written between `from` and `to`, in simulation
 * periods.  Books still on the market in period `to` aren't included (because they might stay on
 * the market).
 */
double book_market_periods(const state::Storage &cs, eris::eris_time_t from, eris::eris_time_t to);

/** Returns the average first-sale-period price of books written in the given period range.  (This
 * is absolute price, not price less marginal cost).
 */
double book_p0(const state::Storage &cs, eris::eris_time_t from, eris::eris_time_t to);

/** Returns the average second-period price of books written in the given period range.  (This is
 * absolute price, not price less marginal cost).  Note that books written in period `to` are not
 * considered (because their 2nd-period price occurs in `to+1`), but books written in period
 * `from-1` are.
 *
 * If there are no suitable books that were on the market for 2+ periods at all, NaN is returned.
 */
double book_p1(const state::Storage &cs, eris::eris_time_t from, eris::eris_time_t to);

/** Returns the average third-period price of books.  (This is absolute price, not price less
 * marginal cost).  Note that books written in period `to` and `to-1` are not included, since their
 * third-period occurs later than `to`, but books written in periods `from-1` and `from-2` are
 * considered.
 *
 * If there are no suitable books that were on the market for 3+ periods at all, NaN is returned.
 */
double book_p2(const state::Storage &cs, eris::eris_time_t from, eris::eris_time_t to);

/** Average copies sold per book.  All books on the market in the given range are included.  The
 * average is per book seen in the period, not per simulation period.
 */
double book_sales(const state::Storage &cs, eris::eris_time_t from, eris::eris_time_t to);

/** Average per-book revenue over the given period.  Books written before `from` are included (if
 * revenue in incurred in [from,to]; for books who continue selling after `to`, only the revenue up
 * to `to` is included.  All books on the market in the given period range (even if sales are 0) are
 * included.
 *
 * The average is calculated based on the number of books seen, not the number of simulation
 * periods.
 */
double book_revenue(const state::Storage &cs, eris::eris_time_t from, eris::eris_time_t to);

/** Average gross margin (i.e. P-MC) of book over the given period.  Books written before `from` are
 * included (if revenue in incurred in [from,to]; for books who continue selling after `to`, only
 * the revenue up to `to` is included.  All books on the market (even if sales are 0) are included.
 *
 * The average is calculated based on the number of books seen, not the number of simulation
 * periods.
 */
double book_gross_margin(const state::Storage &cs, eris::eris_time_t from, eris::eris_time_t to);

/** Average net profit (i.e. profit minus writing cost and keep-on-market costs) of a book.  Only
 * costs incurred during the period are included.  In particular, this means profits from
 * pre-`from`-written books are included but (some) fixed costs are not, and some near-`to` books
 * will have fixed costs but may omit some earned profits.
 *
 * The average is calculated based on the number of books seen, not the number of simulation
 * periods.
 */
double book_profit(const state::Storage &cs, eris::eris_time_t from, eris::eris_time_t to);

/** Average quality of books written during the period range.
 */
double book_quality(const state::Storage &cs, eris::eris_time_t from, eris::eris_time_t to);

/** Average number of books written per period.
 */
double books_written(const state::Storage &cs, eris::eris_time_t from, eris::eris_time_t to);

/** Average number of books purchased per period (aggregate, not per-reader)
 */
double books_bought(const state::Storage &cs, eris::eris_time_t from, eris::eris_time_t to);

/** Average number of books pirated per period (aggreate, not per-reader).  Note that pirated copies
 * of books that left the market before `from` are still included.
 */
double books_pirated(const state::Storage &cs, eris::eris_time_t from, eris::eris_time_t to);

/** Average number of books obtained from the public provider.
 */
double books_public_copies(const state::Storage &cs, eris::eris_time_t from, eris::eris_time_t to);

/** Returns a vector of all supported initial datum values. */
std::vector<initial_datum> initial_data_fields();

/** Returns a vector of all supported calculated datum values. */
std::vector<datum> data_fields();

/** Takes a string, manipulates it into our simplified CSV-suitable value by removing any newlines,
 * commas, and quotation marks.
 */
std::string csv_fix(std::string val);

}}
