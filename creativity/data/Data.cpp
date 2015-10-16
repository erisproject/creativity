#include "creativity/data/Data.hpp"
#include "creativity/state/State.hpp"
#include "creativity/Creativity.hpp"
#include "creativity/Reader.hpp"
#include "creativity/state/Storage.hpp"
#include <boost/algorithm/string/erase.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/iterator/iterator_traits.hpp>
#include <limits>
#include <memory>
#include <stdexcept>
#include <map>
#include <unordered_set>
#include <utility>

using namespace creativity::state;
using namespace eris;

namespace creativity { namespace data {

/** Calculates the average (private) market life of books written between `from` and `to`, in
 * simulation periods.  Books still on the market in period `to` aren't included (because they might
 * stay on the market).
 */
double book_market_periods(const Storage &cs, eris_time_t from, eris_time_t to) {
    if (from > to) throw std::logic_error("from > to");
    unsigned long total = 0;
    unsigned int count = 0;
    auto csto = cs[to];
    for (auto &bp : csto->books) {
        auto &b = bp.second;
        if (b.created >= from and not b.market_private) {
            total += b.lifetime_private;
            count++;
        }
    }

    return total / (double) count;
}

/** Returns the average first-sale-period price of books written in the given period range.  (This
 * is absolute price, not price less marginal cost).
 */
double book_p0(const Storage &cs, eris_time_t from, eris_time_t to) {
    if (from > to) throw std::logic_error("from > to");
    double p_total = 0;
    unsigned int count = 0;
    for (eris_time_t t = std::max<eris_time_t>(1, from); t <= to; t++) {
        auto cst = cs[t];
        for (auto &bp : cst->books) {
            auto &b = bp.second;
            if (b.created == t and b.market_private) {
                p_total += b.price;
                count++;
            }
        }
    }
    if (count == 0) return std::numeric_limits<double>::quiet_NaN();
    return p_total / count;
}

/** Returns the average second-period price of books written in the given period range.  (This is
 * absolute price, not price less marginal cost).  Note that books written in period `to` are not
 * considered (because their 2nd-period price occurs in `to+1`), but books written in period
 * `from-1` are.
 *
 * If there are no suitable books that were on the market for 2+ periods at all, NaN is returned.
 */
double book_p1(const Storage &cs, eris_time_t from, eris_time_t to) {
    if (from > to) throw std::logic_error("from > to");
    double p_total = 0;
    unsigned int count = 0;
    for (eris_time_t t = std::max<eris_time_t>(2, from); t <= to; t++) {
        auto cst = cs[t];
        for (auto &bp : cst->books) {
            auto &b = bp.second;
            if (b.created == t-1 and b.market_private) {
                p_total += b.price;
                count++;
            }
        }
    }
    if (count == 0) return std::numeric_limits<double>::quiet_NaN();
    return p_total / count;
}

/** Returns the average third-period price of books.  (This is absolute price, not price less
 * marginal cost).  Note that books written in period `to` and `to-1` are not included, since their
 * third-period occurs later than `to`, but books written in periods `from-1` and `from-2` are
 * considered.
 *
 * If there are no suitable books that were on the market for 3+ periods at all, NaN is returned.
 */
double book_p2(const Storage &cs, eris_time_t from, eris_time_t to) {
    if (from > to) throw std::logic_error("from > to");
    double p_total = 0;
    unsigned int count = 0;
    for (eris_time_t t = std::max<eris_time_t>(3, from); t <= to; t++) {
        auto cst = cs[t];
        for (auto &bp : cst->books) {
            auto &b = bp.second;
            if (b.created == t-2 and b.market_private) {
                p_total += b.price;
                count++;
            }
        }
    }
    if (count == 0) return std::numeric_limits<double>::quiet_NaN();
    return p_total / count;
}

/** Average private copies sold per book.  All books on the private market in the given range are
 * included.  The average is per book seen in the period, not per simulation period.
 */
double book_sales(const Storage &cs, eris_time_t from, eris_time_t to) {
    if (from > to) throw std::logic_error("from > to");
    long sales_total = 0;
    std::unordered_set<eris_id_t> seen;
    for (eris_time_t t = from; t <= to; t++) {
        auto cst = cs[t];
        for (auto &bp : cst->books) {
            auto &b = bp.second;
            if (b.market_private) {
                sales_total += b.sales;
                seen.insert(b.id);
            }
        }
    }
    return sales_total / (double) seen.size();
}

/** Average per-book revenue over the given period.  Books written before `from` are included (if
 * revenue in incurred in [from,to]; for books who continue selling after `to`, only the revenue up
 * to `to` is included.  All books on the market in the given period range (even if sales are 0) are
 * included.
 *
 * The average is calculated based on the number of books seen, not the number of simulation
 * periods.
 */
double book_revenue(const Storage &cs, eris_time_t from, eris_time_t to) {
    if (from > to) throw std::logic_error("from > to");
    double revenue_total = 0;
    std::unordered_set<eris_id_t> seen;
    for (eris_time_t t = from; t <= to; t++) {
        auto cst = cs[t];
        for (auto &bp : cst->books) {
            auto &b = bp.second;
            if (b.market_private) {
                revenue_total += b.revenue;
                seen.insert(b.id);
            }
        }
    }

    return revenue_total / seen.size();
}

/** Average gross margin (i.e. P-MC) of book over the given period.  Books written before `from` are
 * included (if revenue in incurred in [from,to]; for books who continue selling after `to`, only
 * the revenue up to `to` is included.  All books on the market (even if sales are 0) are included.
 *
 * The average is calculated based on the number of books seen, not the number of simulation
 * periods.
 */
double book_gross_margin(const Storage &cs, eris_time_t from, eris_time_t to) {
    if (from > to) throw std::logic_error("from > to");
    double margin_total = 0;
    std::unordered_set<eris_id_t> seen;
    for (eris_time_t t = from; t <= to; t++) {
        auto cst = cs[t];
        for (const auto &bp : cst->books) {
            auto &b = bp.second;
            if (b.market_private) {
                margin_total += b.revenue -  cs.settings.cost_unit * b.sales;
                seen.insert(b.id);
            }
        }
    }

    return margin_total / seen.size();
}

/** Average net profit (i.e. profit minus writing cost and keep-on-market costs) of a book.  Only
 * costs incurred during the period are included.  In particular, this means profits from
 * pre-`from`-written books are included but (some) fixed costs are not, and some near-`to` books
 * will have fixed costs but may omit some earned profits.
 *
 * The average is calculated based on the number of books seen, not the number of simulation
 * periods.
 */
double book_profit(const Storage &cs, eris_time_t from, eris_time_t to) {
    if (from > to) throw std::logic_error("from > to");
    double profit_total = 0;
    std::unordered_set<eris_id_t> seen;
    for (eris_time_t t = from; t <= to; t++) {
        auto period = cs[t];
        for (const auto &bp : period->books) {
            auto &b = bp.second;
            if (b.market_private) {
                seen.insert(b.id);
                profit_total += b.revenue - cs.settings.cost_unit * b.sales - cs.settings.cost_market;
                if (b.created == t) {
                    // In the creation period, subtract the effort that had to be expended to write
                    // the book
                    auto &r = period->readers.at(b.author);
                    profit_total -= Reader::creationEffort(r.creation_shape, r.creation_scale, b.quality);
                }
            }
        }
    }

    return profit_total / seen.size();
}

struct quantile_cache {
    eris_time_t from = 0, to = 0;
    const Storage *cs = nullptr;
    std::vector<double> values;
    double mean = std::numeric_limits<double>::quiet_NaN();
};
thread_local std::unordered_map<std::string, quantile_cache> q_cache;

// Loads a quantile_cache with values calculated from books created in the given period range
const quantile_cache& load_quantile_cache(const std::string &variable, const Storage &cs, eris_time_t from, eris_time_t to) {
    if (from > to) throw std::logic_error("from > to");
    auto &cache = q_cache[variable];
    if (cache.cs == &cs and cache.from == from and cache.to == to) {
        return cache; // Cache is good
    }

    // We either loop through books, adding a data point for each book written in the period, or we
    // loop through readers, calculating an average value for each reader through the period, then
    // calculate the sd and quantiles from the per-reader temporal means.
    bool books;
    std::function<double(const BookState &b)> book_val;
    std::function<double(const ReaderState &r)> reader_val;
    if (variable == "book_quality") {
        books = true;
        book_val = [](const BookState &b) { return b.quality; };
    }
    else if (variable == "book_author_scale") {
        books = true;
        book_val = [&cs](const BookState &b) { return cs[b.created]->readers.at(b.author).creation_scale; };
    }
    else if (variable == "book_author_effort") {
        books = true;
        book_val = [&cs](const BookState &b) {
            auto &a = cs[b.created]->readers.at(b.author);
            return Reader::creationEffort(a.creation_shape, a.creation_scale, b.quality);
        };
    }
    else if (variable == "net_u") {
        books = false;
        reader_val = [&cs](const ReaderState &r) { return r.u - cs.settings.income; };
    }
    else {
        throw std::logic_error("Unknown quantile variable: " + variable);
    }

    cache.values.clear();
    cache.cs = &cs;
    cache.from = from;
    cache.to = to;
    cache.mean = 0;

    if (books) {
        for (eris_time_t t = from; t <= to; t++) {
            for (auto &bp : cs[t]->books) {
                auto &b = bp.second;
                if (b.created == t) {
                    cache.values.push_back(book_val(b));
                    cache.mean += cache.values.back();
                }
            }
        }
    }
    else {
        std::map<eris_id_t, std::pair<double, unsigned>> reader_sum;
        for (eris_time_t t = from; t <= to; t++) {
            for (auto &rp : cs[t]->readers) {
                auto &rm = reader_sum[rp.second.id];
                rm.first += reader_val(rp.second);
                rm.second++;
            }
        }
        for (auto &id_val : reader_sum) {
            cache.values.push_back(id_val.second.first / id_val.second.second);
            cache.mean += cache.values.back();
        }
    }
    if (cache.values.empty()) cache.mean = std::numeric_limits<double>::quiet_NaN();
    else cache.mean /= cache.values.size();

    std::sort(cache.values.begin(), cache.values.end());

    return cache;
}

double quantile(const std::vector<double> &vals, double prob) {
    return quantile(Eigen::Map<const Eigen::VectorXd>(vals.data(), vals.size()), prob);
}

double quantile(const Eigen::Ref<const Eigen::VectorXd> &vals, double prob) {
    if (prob < 0 or prob > 1) throw std::logic_error("Requested quantile probability is invalid");
    if (vals.size() == 0) return std::numeric_limits<double>::quiet_NaN();
    double index = prob * (vals.size()-1);
    unsigned below = std::floor(index), above = std::ceil(index);
    if (below == above) return vals[above];
    return (above - index) * vals[below] + (index - below) * vals[above];
}

double variance(const std::vector<double> &vals, double mean) {
    if (vals.size() < 2) return std::numeric_limits<double>::quiet_NaN();
    if (std::isnan(mean)) {
        mean = 0;
        for (const auto &v : vals) mean += v;
        mean /= vals.size();
    }
    // Keep a sum of the demeaned values: mathematically they should equal 0, but numerical error
    // can creep in, and so this is essentially collecting the numerical error so that it can be
    // corrected for.
    double demeaned_squares = 0, demeaned_sum = 0;
    for (const auto &v : vals) {
        double d = v - mean;
        demeaned_squares += d*d;
        demeaned_sum += d;
    }
    return (demeaned_squares - demeaned_sum*demeaned_sum / vals.size()) / (vals.size()-1);
}

#define MEAN_FN(variable) double variable(const Storage &cs, eris_time_t from, eris_time_t to) { \
    return load_quantile_cache(#variable, cs, from, to).mean; \
}
#define SD_FN(variable) double variable##_sd(const Storage &cs, eris_time_t from, eris_time_t to) { \
    const auto &cache = load_quantile_cache(#variable, cs, from, to); \
    if (cache.values.size() < 2) return std::numeric_limits<double>::quiet_NaN(); \
    return std::sqrt(variance(cache.values, cache.mean)); \
}
#define QUANTILE_FN(variable, suffix, prob) double variable##_##suffix(const state::Storage &cs, eris_time_t from, eris_time_t to) { \
    return quantile(load_quantile_cache(#variable, cs, from, to).values, prob); \
}
#define DIST_FNS(variable) \
MEAN_FN(variable) \
SD_FN(variable) \
QUANTILE_FN(variable, min, 0) \
QUANTILE_FN(variable, 5th, 0.05) \
QUANTILE_FN(variable, 10th, 0.1) \
QUANTILE_FN(variable, 25th, 0.25) \
QUANTILE_FN(variable, median, 0.5) \
QUANTILE_FN(variable, 75th, 0.75) \
QUANTILE_FN(variable, 90th, 0.9) \
QUANTILE_FN(variable, 95th, 0.95) \
QUANTILE_FN(variable, max, 1)

DIST_FNS(net_u)

DIST_FNS(book_quality)

DIST_FNS(book_author_scale)

DIST_FNS(book_author_effort)

#undef MEAN_FN
#undef SD_FN
#undef QUANTILE_FN
#undef DIST_FNS

/** Average number of books written per period.
 */
double books_written(const Storage &cs, eris_time_t from, eris_time_t to) {
    if (from > to) throw std::logic_error("from > to");
    unsigned int count = 0;
    for (eris_time_t t = from; t <= to; t++) {
        auto cst = cs[t];
        for (auto &bp : cst->books) {
            auto &b = bp.second;
            if (b.created == t) {
                count++;
            }
        }
    }

    return count / (double) (to-from+1);
}

/** Average number of books purchased privately per period (aggregate, not per-reader)
 */
double books_bought(const Storage &cs, eris_time_t from, eris_time_t to) {
    if (from > to) throw std::logic_error("from > to");
    unsigned long count = 0;
    for (eris_time_t t = from; t <= to; t++) {
        auto cst = cs[t];
        for (auto &bp : cst->books) {
            auto &b = bp.second;

            if (b.market_private)
                count += b.sales;
        }
    }

    return count / (double) (to-from+1);
}

/** Average number of books pirated per period (aggreate, not per-reader).  Note that pirated copies
 * of books that left the market before `from` are still included.
 */
double books_pirated(const Storage &cs, eris_time_t from, eris_time_t to) {
    if (from > to) throw std::logic_error("from > to");
    unsigned long count = 0;
    for (eris_time_t t = from; t <= to; t++) {
        auto cst = cs[t];
        for (auto &bp : cst->books) {
            auto &b = bp.second;

            count += b.pirated;
        }
    }

    return count / (double) (to-from+1);
}

/** Average number of public copies of books provided per period (aggreate, not per-reader).
 */
double books_public_copies(const Storage &cs, eris_time_t from, eris_time_t to) {
    if (from > to) throw std::logic_error("from > to");
    unsigned long count = 0;
    for (eris_time_t t = from; t <= to; t++) {
        auto cst = cs[t];
        for (auto &bp : cst->books) {
            auto &b = bp.second;

            if (b.market_public())
                count += b.sales;
        }
    }

    return count / (double) (to-from+1);
}

double reader_market_spending(const Storage &cs, eris_time_t from, eris_time_t to) {
    if (from > to) throw std::logic_error("from > to");
    double avg_spending = 0;
    for (eris_time_t t = from; t <= to; t++) {
        auto cst = cs[t];
        double period_spending = 0;
        auto &books = cst->books;
        for (auto &rdr : cst->readers) {
            for (auto &bkid : rdr.second.new_books) {
                if (rdr.second.library.at(bkid).purchased_market()) {
                    period_spending += books.at(bkid).price;
                }
            }
        }
        avg_spending += period_spending / cst->readers.size();
    }
    return avg_spending / (to-from+1);
}

double reader_piracy_spending(const Storage &cs, eris_time_t from, eris_time_t to) {
    if (from > to) throw std::logic_error("from > to");
    double avg_spending = 0;
    for (eris_time_t t = from; t <= to; t++) {
        auto cst = cs[t];
        double period_spending = 0;
        for (auto &rdr : cst->readers) {
            for (auto &bkid : rdr.second.new_books) {
                if (rdr.second.library.at(bkid).pirated()) {
                    period_spending += cs.settings.cost_piracy;
                }
            }
        }
        avg_spending += period_spending / cst->readers.size();
    }
    return avg_spending / (to-from+1);
}

double reader_public_spending(const Storage &cs, eris_time_t from, eris_time_t to) {
    if (from > to) throw std::logic_error("from > to");
    double avg_spending = 0;
    for (eris_time_t t = from; t <= to; t++) {
        auto cst = cs[t];
        double period_spending = 0;
        auto &books = cst->books;
        for (auto &rdr : cst->readers) {
            for (auto &bkid : rdr.second.new_books) {
                if (rdr.second.library.at(bkid).purchased_public()) {
                    period_spending += books.at(bkid).price;
                }
            }
        }
        avg_spending += period_spending / cst->readers.size();
    }
    return avg_spending / (to-from+1);
}

double reader_spending(const Storage &cs, eris_time_t from, eris_time_t to) {
    if (from > to) throw std::logic_error("from > to");
    double avg_spending = 0;
    for (eris_time_t t = from; t <= to; t++) {
        auto cst = cs[t];
        double period_spending = 0;
        auto &books = cst->books;
        for (auto &rdr : cst->readers) {
            for (auto &bkid : rdr.second.new_books) {
                auto &bookcopy = rdr.second.library.at(bkid);
                period_spending +=
                    (bookcopy.purchased_market() or bookcopy.purchased_public())
                    ? books.at(bkid).price
                    : cs.settings.cost_piracy;
            }
        }
        avg_spending += period_spending / cst->readers.size();
    }
    return avg_spending / (to-from+1);
}


std::vector<initial_datum> initial_data_fields() {
    std::vector<initial_datum> initial_data;
#define ADD_SETTING(S) initial_data.emplace_back(#S, [](const CreativitySettings &cs) { return cs.S; })
    ADD_SETTING(readers);
    ADD_SETTING(dimensions);
    ADD_SETTING(boundary);
    initial_data.emplace_back("density", [](const CreativitySettings &cs) { return Creativity::densityFromBoundary(cs.readers, cs.dimensions, cs.boundary); });
    ADD_SETTING(book_distance_mean);
    ADD_SETTING(book_quality_sd);
    ADD_SETTING(reader_step_mean);
    ADD_SETTING(reader_creation_shape);
    ADD_SETTING(reader_creation_scale_min);
    ADD_SETTING(reader_creation_scale_range);
    ADD_SETTING(cost_market);
    ADD_SETTING(cost_unit);
    ADD_SETTING(cost_piracy);
    ADD_SETTING(creation_time);
    ADD_SETTING(creation_fixed);
    ADD_SETTING(income);
    ADD_SETTING(piracy_begins);
    ADD_SETTING(piracy_link_proportion);
    ADD_SETTING(public_sharing_begins);
    ADD_SETTING(public_sharing_tax);
    ADD_SETTING(prior_scale);
    ADD_SETTING(prior_scale_piracy);
    ADD_SETTING(prior_scale_burnin);
    ADD_SETTING(burnin_periods);
    ADD_SETTING(prediction_draws);
    ADD_SETTING(initial.prob_write);
    ADD_SETTING(initial.l_min);
    ADD_SETTING(initial.l_range);
    ADD_SETTING(initial.p_min);
    ADD_SETTING(initial.p_range);
    ADD_SETTING(initial.prob_keep);
    ADD_SETTING(initial.keep_price);
    ADD_SETTING(initial.belief_threshold);
#undef ADD_SETTING

    return initial_data;
}

std::vector<datum> data_fields() {
    std::vector<datum> data;

#define ADD_DATUM(D, ...) data.emplace_back(#D, D, ##__VA_ARGS__)
#define ADD_DIST_DATA(variable) \
    ADD_DATUM(variable); \
    ADD_DATUM(variable##_sd); \
    ADD_DATUM(variable##_min); \
    ADD_DATUM(variable##_5th); \
    ADD_DATUM(variable##_10th); \
    ADD_DATUM(variable##_25th); \
    ADD_DATUM(variable##_median); \
    ADD_DATUM(variable##_75th); \
    ADD_DATUM(variable##_90th); \
    ADD_DATUM(variable##_95th); \
    ADD_DATUM(variable##_max);
    ADD_DIST_DATA(net_u);
    ADD_DATUM(book_market_periods);
    ADD_DATUM(book_p0);
    ADD_DATUM(book_p1);
    ADD_DATUM(book_p2);
    ADD_DATUM(book_sales);
    ADD_DATUM(book_revenue);
    ADD_DATUM(book_gross_margin);
    ADD_DATUM(book_profit);
    ADD_DIST_DATA(book_quality);
    ADD_DATUM(books_written);
    ADD_DATUM(books_bought);
    ADD_DATUM(books_pirated, false, true, true);
    ADD_DATUM(books_public_copies, false, false, true);
    ADD_DATUM(reader_spending);
    ADD_DATUM(reader_market_spending);
    ADD_DATUM(reader_piracy_spending, false, true, true);
    ADD_DATUM(reader_public_spending, false, false, true);
    ADD_DIST_DATA(book_author_scale);
    ADD_DIST_DATA(book_author_effort);
#undef ADD_DATUM

    return data;
}

std::string csv_fix(std::string val) {
    boost::erase_all(val, "\n"); // remove newlines
    boost::erase_all(val, "\""); // remove "s
    boost::erase_all(val, ","); // remove ,s
    return val;
}

}}
