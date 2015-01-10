#include "creativity/data/Data.hpp"
#include "creativity/Creativity.hpp"
#include <boost/algorithm/string.hpp>

using namespace creativity::state;
using namespace eris;

namespace creativity { namespace data {

/** Calculates the average net utility over the given period.  Net utility is a reader's utility
 * minus the reader's income (since without any simulation activity, the reader's utility is
 * quasilinear, thus the reader receives exactly his income as utility).
 */
double net_u(const Storage &cs, eris_time_t from, eris_time_t to) {
    if (from > to) throw std::logic_error("from > to");
    double net_u_total = 0;
    for (eris_time_t t = from; t <= to; t++) {
        auto cst = cs[t];
        for (auto &r : cst->readers) {
            net_u_total += r.second.u - r.second.income;
        }
    }
    return net_u_total / (cs[from]->readers.size() * (from-to+1));
}

/** Calculates the average market life of books written between `from` and `to`, in simulation
 * periods.  Books still on the market in period `to` aren't included (because they might stay on
 * the market).
 */
double book_market_periods(const Storage &cs, eris_time_t from, eris_time_t to) {
    if (from > to) throw std::logic_error("from > to");
    unsigned long total = 0;
    unsigned int count = 0;
    auto csto = cs[to];
    for (auto &bp : csto->books) {
        auto &b = bp.second;
        if (b.created >= from and not b.market()) {
            total += b.lifetime;
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
            if (b.created == t and b.market()) {
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
            if (b.created == t-1 and b.market()) {
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
            if (b.created == t-2 and b.market()) {
                p_total += b.price;
                count++;
            }
        }
    }
    if (count == 0) return std::numeric_limits<double>::quiet_NaN();
    return p_total / count;
}

/** Average copies sold per book.  All books on the market in the given range are included.  The
 * average is per book seen in the period, not per simulation period.
 */
double book_sales(const Storage &cs, eris_time_t from, eris_time_t to) {
    if (from > to) throw std::logic_error("from > to");
    long sales_total = 0;
    std::unordered_set<eris_id_t> seen;
    for (eris_time_t t = from; t <= to; t++) {
        auto cst = cs[t];
        for (auto &bp : cst->books) {
            auto &b = bp.second;
            if (b.market()) {
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
            if (b.market()) {
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
            if (b.market()) {
                auto &r = cst->readers.at(b.author);
                margin_total += b.revenue - r.cost_unit * b.sales;
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
            if (b.market()) {
                seen.insert(b.id);
                auto &r = period->readers.at(b.author);
                profit_total += b.revenue - r.cost_unit * b.sales - r.cost_fixed;
                if (b.created == t) {
                    // In the creation period, subtract the effort that had to be expended to write
                    // the book
                    profit_total -= Reader::creationEffort(r.creation_shape, r.creation_scale, b.quality);
                }
            }
        }
    }

    return profit_total / seen.size();
}

/** Average quality of books written during the period range.
 */
double book_quality(const Storage &cs, eris_time_t from, eris_time_t to) {
    if (from > to) throw std::logic_error("from > to");
    double q_total = 0;
    unsigned int count = 0;
    for (eris_time_t t = from; t <= to; t++) {
        auto cst = cs[t];
        for (auto &bp : cst->books) {
            auto &b = bp.second;
            if (b.created == t) {
                q_total += b.quality;
                count++;
            }
        }
    }

    return q_total / count;
}

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

/** Average number of books purchased per period (aggregate, not per-reader)
 */
double books_bought(const Storage &cs, eris_time_t from, eris_time_t to) {
    if (from > to) throw std::logic_error("from > to");
    unsigned long count = 0;
    for (eris_time_t t = from; t <= to; t++) {
        auto cst = cs[t];
        for (auto &bp : cst->books) {
            auto &b = bp.second;

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

std::vector<initial_datum> initial_data_fields() {
    std::vector<initial_datum> initial_data;
#define ADD_SETTING(S) initial_data.emplace_back(#S, [](const CreativitySettings &cs) { return cs.S; })
    ADD_SETTING(readers);
    ADD_SETTING(dimensions);
    ADD_SETTING(boundary);
    initial_data.emplace_back("density", [](const CreativitySettings &cs) { return Creativity::densityFromBoundary(cs.readers, cs.dimensions, cs.boundary); });
    ADD_SETTING(book_distance_sd);
    ADD_SETTING(book_quality_sd);
    ADD_SETTING(reader_step_sd);
    ADD_SETTING(reader_creation_shape);
    ADD_SETTING(reader_creation_scale_min);
    ADD_SETTING(reader_creation_scale_max);
    ADD_SETTING(cost_fixed);
    ADD_SETTING(cost_unit);
    ADD_SETTING(cost_piracy);
    ADD_SETTING(income);
    ADD_SETTING(piracy_begins);
    ADD_SETTING(piracy_link_proportion);
    ADD_SETTING(prior_weight);
    ADD_SETTING(prior_weight_piracy);
    ADD_SETTING(initial.prob_write);
    ADD_SETTING(initial.q_min);
    ADD_SETTING(initial.q_max);
    ADD_SETTING(initial.p_min);
    ADD_SETTING(initial.p_max);
    ADD_SETTING(initial.prob_keep);
    ADD_SETTING(initial.keep_price);
    ADD_SETTING(initial.belief_threshold);
#undef ADD_SETTING

    return initial_data;
}

std::vector<datum> data_fields() {
    std::vector<datum> data;

#define ADD_DATUM(D) data.emplace_back(#D, D)
    ADD_DATUM(net_u);
    ADD_DATUM(book_market_periods);
    ADD_DATUM(book_p0);
    ADD_DATUM(book_p1);
    ADD_DATUM(book_p2);
    ADD_DATUM(book_sales);
    ADD_DATUM(book_revenue);
    ADD_DATUM(book_gross_margin);
    ADD_DATUM(book_profit);
    ADD_DATUM(book_quality);
    ADD_DATUM(books_written);
    ADD_DATUM(books_bought);
    data.emplace_back("books_pirated", books_pirated, true);
#undef ADD_DATUM

    return data;
}

std::string csv_escape(std::string val) {
    // First, double-up any quotation markseliminate newlines.  It seems pretty unlikely that the input could contain them, but
    // just in case:
    boost::erase_all(val, "\n"); // remove newlines
    boost::replace_all(val, "\"", "\"\""); // double-up quotes
    // If there's a "", comma, or semicolon then quote the whole thing (the semicolon isn't really
    // needed, but might help not confusing a CSV parser).
    if (boost::algorithm::contains(val, "\"\"") or boost::algorithm::contains(val, ",") or boost::algorithm::contains(val, ";"))
        return "\"" + val + "\"";
    else
        return val;
}

}}
