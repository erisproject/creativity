#pragma once
#include <string>
#include <functional>
#include <vector>
#include <algorithm>
#include <type_traits>
#include <eris/types.hpp>
#include <Eigen/Core>

/** \file
 * Various routines for generating data from stored creativity results.
 */

namespace creativity { struct CreativitySettings; }
namespace creativity { namespace state { class Storage; } }

namespace creativity {
/// Namespace for classes related to generating and processing simulation data.
namespace data {

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
    initial_datum(std::string name, std::function<IntType(const CreativitySettings&)> &&calc_int)
        : name(std::move(name)), calc_int(std::move(calc_int)) {}
    /// Constructs a simulation parameter datum that returns a double value.
    initial_datum(std::string name, std::function<double(const CreativitySettings&)> &&calc_double)
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
        bool policy; ///< True if this field applies to policy periods
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
     * \param policy specifies whether this data field applies to policy periods (default true)
     */
    datum(std::string &&name, decltype(calculate) calc, bool pre = true, bool piracy = true, bool policy = true)
        : name(std::move(name)), calculate(std::move(calc)), applies_to{pre, piracy, policy} {}
};

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
 * revenue is received in [from,to]; for books who continue selling after `to`, only the revenue up
 * to `to` is included.  All books created in the given period range (even if revenue is 0) are
 * included.
 *
 * Prize revenue (i.e. public market payouts) is counted as revenue.
 *
 * The average is calculated based on the number of books included, not the number of simulation
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

/** Average profit (i.e. profit minus writing cost and keep-on-market costs) of a book.  Only costs
 * incurred during the period are included.  In particular, this means profits from
 * pre-`from`-written books are included but (some) fixed costs are not, and some near-`to` books
 * will have fixed costs but may omit some earned profits.
 *
 * Prize revenue (i.e. public market payours) are included.
 *
 * The average is calculated based on the number of books seen, not the number of simulation
 * periods.
 */
double book_profit(const state::Storage &cs, eris::eris_time_t from, eris::eris_time_t to);

/// Macro to generate mean, sd, and quantile functions for a given variable
#define DIST_FNS(variable) \
double variable(const state::Storage &cs, eris::eris_time_t from, eris::eris_time_t to); \
double variable##_sd(const state::Storage &cs, eris::eris_time_t from, eris::eris_time_t to); \
double variable##_min(const state::Storage &cs, eris::eris_time_t from, eris::eris_time_t to); \
double variable##_5th(const state::Storage &cs, eris::eris_time_t from, eris::eris_time_t to); \
double variable##_10th(const state::Storage &cs, eris::eris_time_t from, eris::eris_time_t to); \
double variable##_25th(const state::Storage &cs, eris::eris_time_t from, eris::eris_time_t to); \
double variable##_median(const state::Storage &cs, eris::eris_time_t from, eris::eris_time_t to); \
double variable##_75th(const state::Storage &cs, eris::eris_time_t from, eris::eris_time_t to); \
double variable##_90th(const state::Storage &cs, eris::eris_time_t from, eris::eris_time_t to); \
double variable##_95th(const state::Storage &cs, eris::eris_time_t from, eris::eris_time_t to); \
double variable##_max(const state::Storage &cs, eris::eris_time_t from, eris::eris_time_t to);
///@{
/** Distribution values (mean, standard deviation, min, max, quantiles) for quality of books written
 * during the period range.
 */
DIST_FNS(book_quality)
///@}

///@{
/** Distribution values (mean, standard deviation, min, max, quantiles) for author creation scale of
 * the author of each book created during the period range.
 */
DIST_FNS(book_author_level)
///@}

///@{
/** Distribution values (mean, standard deviation, min, max, quantiles) for author effort level of
 * each book created during the period range.
 */
DIST_FNS(book_author_effort)
///@}

///@{
/** Distribution of reader net utility values over the period range.  This is net of the utility
 * readers would have if there were no book activity at all, which is simply the reader's exogenous
 * income level.
 */
DIST_FNS(net_u)
///@}

#undef DIST_FNS

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

/** Average per-reader private market spending on books, averaged over the given period.  Only
 * spending on books bought from the author is included, specifically piracy cost and public sharing
 * costs are not included.
 */
double reader_market_spending(const state::Storage &cs, eris::eris_time_t from, eris::eris_time_t to);

/** Average per-reader piracy cost expenditure for obtaining pirated copies of books, averaged over
 * the given period.
 */
double reader_piracy_spending(const state::Storage &cs, eris::eris_time_t from, eris::eris_time_t to);

/** Average per-reader public sharing expenditure for obtaining public provider copies of books,
 * averaged over the given period.
 *
 * This does *not* include the lump-sum tax paid by all readers, only the per-copy price (= marginal
 * cost) of a book.
 */
double reader_public_spending(const state::Storage &cs, eris::eris_time_t from, eris::eris_time_t to);

/** Average per-reader expenditure for obtaining books from any source.  This is the sum of
 * book_market_spending plus book_piracy_spending (when piracy is available) plus
 * book_public_spending (when public sharing is available).
 *
 * This does *not* include spending incurred by authors to create books or sell copies of books, nor
 * does it include lump sum tax amounts.
 */
double reader_spending(const state::Storage &cs, eris::eris_time_t from, eris::eris_time_t to);


/** Returns a vector of all supported initial datum values. */
std::vector<initial_datum> initial_data_fields();

/** Returns a vector of all supported calculated datum values. */
std::vector<datum> data_fields();

/** Takes a string, manipulates it into our simplified CSV-suitable value by replacing any special
 * characters (such as commas, newline, etc.) with underscores.
 */
std::string csv_fix(std::string val);

}}
