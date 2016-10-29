#pragma once
#include <creativity/data/Treatment.hpp>
#include <memory>

namespace creativity { namespace data {

/** Extension to Treatment that applies a filter to an existing Treatment variable, resulting in a
 * new Treatment object that only contains the row sets for matched data observations.
 *
 * For example, the following:
 *
 *     Treatment data("some.csv");
 *     TreatmentFilter data_filtered(data, [](const TreatmentFilter::Properties &p) {
 *         return p.piracy and p.piracy->value("books_written") >= 0.2;
 *     });
 *
 * constructs a filter that only includes data observations that have piracy data with at a piracy
 * books_written value of at least 0.2.
 *
 * The filter is applied on a simulation-by-simulation row, which typically corresponds to multiple
 * treatment rows at once.
 *
 * The filter can also be applied to select only particular types of rows: for example, the
 * following constructs a filter that only includes "pre" and "piracy" rows, but omits short-run
 * piracy, and public rows:
 *
 *     Treatment data("some.csv");
 *     TreatmentFilter pre_and_piracy(
 *         data,
 *         [](const TreatmentFilter::Properties&) { return true; } // No per-simulation filtering
 *         [](bool pre, bool piracy, bool policy, bool short_run) {
 *             return pre or (piracy and not short_run);
 *         }
 *     );
 */
class TreatmentFilter : public Treatment {
    public:
        /// Provides a simple interface to accessing the value of a field for a single observation.
        class StageProperties final {
            public:
                /** Returns the value associated with the given stage field name.  The field should
                 * not be prefixed with the "pre.", "piracy.", etc. prefix, but should be prefixed
                 * with "param." to obtain parameter values.
                 *
                 * \sa Treatment::column
                 * \sa Treatment::variable
                 *
                 * \throws std::out_of_range if the given field does not exist.
                 */
                double value(const std::string &field) const;

            private:
                const Treatment &source_;
                const unsigned row_;
                StageProperties(const Treatment &source, unsigned rownum)
                    : source_(source), row_(rownum) {}
                friend class TreatmentFilter;
        };

        /// Data properties for rows filtered by TreatmentFilter()
        struct Properties final {
            public:
                /// A pointer to the Properties value for the pre row associated with this row.
                std::unique_ptr<const StageProperties> pre;
                /** A pointer to the Properties values for the LR piracy row associated with this
                 * row.  Will be a nullptr if there is no such row. Note that whether this is set
                 * depends on whether or not the field is in the source data, *not* whether or not
                 * the stage will be included in the filtered data. */
                std::unique_ptr<const StageProperties> piracy;
                /** A pointer to the Properties values for the LR policy row associated with this
                 * row.  Will be a nullptr if there is no such row. Note that whether this is set
                 * depends on whether or not the field is in the source data, *not* whether or not
                 * the stage will be included in the filtered data. */
                std::unique_ptr<const StageProperties> policy;
                /** A pointer to the Properties values for the SR piracy row associated with this
                 * row.  Will be a nullptr if there is no such row. Note that whether this is set
                 * depends on whether or not the field is in the source data, *not* whether or not
                 * the stage will be included in the filtered data. */
                std::unique_ptr<const StageProperties> piracy_SR;
                /** A pointer to the Properties values for the SR public row associated with this
                 * row.  Will be a nullptr if there is no such row. Note that whether this is set
                 * depends on whether or not the field is in the source data, *not* whether or not
                 * the stage will be included in the filtered data. */
                std::unique_ptr<const StageProperties> policy_SR;
                /// The source data simulation filename.
                std::string source;
        };

        /** Creates a TreatmentFilter object by copying data matching a particular criterion from
         * another TreatmentFilter object.
         *
         * \param source the source Treatment object.
         * \param filter the filter to apply.  This object should return true for any data row sets
         * that should be included in the filtered data.
         * \param stage_filter this filter is called to determine which stages should be considered.
         * It is called multiple times before looping through rows to determine which data source
         * row sets should be passed to `filter`.
         *
         * \throws std::logic_error if the stage filter doesn't admit any stages in `source`
         */
        TreatmentFilter(const Treatment &source,
                std::function<bool(const Properties&)> filter,
                std::function<bool(bool pre, bool piracy, bool policy, bool short_run)> stage_filter = [](bool,bool,bool,bool) { return true; }
                );
};

}}
