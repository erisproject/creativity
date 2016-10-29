#pragma once
#include <creativity/data/CSVParser.hpp>
#include <creativity/data/Variable.hpp>
#include <Eigen/Core>
#include <memory>

namespace creativity { namespace data {

/** This class converts raw simulation data into multiple data rows, with piracy/policy, SR/LR as
 * treatment effects on the base (pre-piracy, pre-policy) row.  It additionally supports filtering
 * the resulting data--for instance, only using observations that still exhibit some writing during
 * piracy.
 */
class Treatment {
    public:
        /// Default constructor
        Treatment() = default;

        /// Takes over the given CSVParser and reads the data from it.
        explicit Treatment(CSVParser &&csv);

        /// Creates a new CSVParser to read the given filename.
        explicit Treatment(const std::string &filename);

        /// Reads the data from the given CSVParser
        void readCSV(CSVParser &&csv);

        /** Requires that the parsed data contain piracy observations.  If the data has already been
         * parsed, this throws an exception immediately if the data does not contain piracy; if not
         * yet parsed, this causes the next readCSV() to throw an exception if piracy data is not
         * present.
         *
         * \param require defaults to true of omitted; can be explicitly specified as false to
         * cancel a previous requirePiracy() call.
         */
        void requirePiracy(bool require = true);

        /// Like requirePiracy(), but for policy data
        void requirePolicy(bool require = true);

        /// Like requirePiracy(), but for pre-piracy data
        void requirePre(bool require = true);

        /** Requires short-run data for long-run piracy/policy data that exists in the source data.
         * In particular, if this option is enabled, the data must contain short-run observations
         * for each category with equivalent long-run observations: that is, if there is long-run
         * piracy data, there must also be short-run piracy data, and likewise for policy data.  If
         * the source data does not contain long-run piracy data, this option will not require
         * short-run piracy data.
         *
         * Like requirePiracy(), this throws immediately if data has already been parsed; if not,
         * the exception will be raised if attempting to read data that doesn't contain the required
         * short-run data.
         */
        void requireSR(bool require = true);

        /// True if this treatment data pre-piracy observations
        const bool& hasPre() const { return has_pre_; }

        /// True if the source data has piracy data
        const bool& hasPiracy() const { return has_piracy_; }

        /// True if the source data has policy data
        const bool& hasPolicy() const { return has_policy_; }

        /// True if the source data has short-run piracy data
        const bool& hasPiracySR() const { return has_piracy_sr_; }

        /// True if the source data has short-run policy data
        const bool& hasPolicySR() const { return has_policy_sr_; }

        /** True if the source data has short-run data for each type of associated long-run data.
         * That is, if the data has piracy data, it must also have short-run piracy data; likewise
         * for policy data.
         */
        bool hasShortrun() const { return (hasPolicySR() or not hasPolicy()) and (hasPiracySR() or not hasPiracy()); }

        /// The number of data rows per simulation input.  This equals 1 plus the number of treatments.
        const unsigned int& rowsPerSimulation() const { return rows_per_sim_; }

        /** The number of simulation files contributing to this data: this is simply an alias for
         * `obj.data().rows() / obj.rowsPerSimulation()`
         */
        unsigned int simulations() const { return data().rows() / rowsPerSimulation(); }

        /** The full set of data.
         *
         * \sa TreatmentFilter for extracting conditional data.
         */
        const Eigen::MatrixXd& data() const { return data_; }

        /** The simulation source values.  Data in row `x` is from source file at index `x /
         * rowsPerSimulation()`.
         */
        const std::vector<std::string>& sourceFiles() const { return source_; }

        /** Return the source file associated with row r.  This is simply an alias for
         * `sourceFiles().at(r / rowsPerSimulation())`.
         */
        const std::string& sourceFile(unsigned r) const { return sourceFiles().at(r / rowsPerSimulation()); }

        /** Accesses the map of field name to column indices.  Data fields (e.g. "books_written")
         * are not prefixed with "pre.", "piracy.", etc.; parameters (e.g. "param.readers") are
         * prefixed; there are also dummies "piracy", "policy", "SR", and "LR".
         */
        const std::unordered_map<std::string, unsigned>& columns() const;

        /** Returns a shared_ptr to a SimpleVariable containing a copy of the requested column.  The
         * column naming values are the same as those accepted by column().
         *
         * Values are internally cached--multiple calls to variable with the same value return the
         * same SimpleVariable.
         *
         * \throws std::out_of_range if the requested variable does not exist.
         */
        std::shared_ptr<const SimpleVariable> variable(const std::string &name) const;

        /// Alias for variable()
        std::shared_ptr<const SimpleVariable> operator[](const std::string &name) const { return variable(name); }

    protected:
        Eigen::MatrixXd data_; ///< The data matrix.
        std::vector<std::string> source_; ///< In-order source values for data.

        bool have_data_{false}, ///< True once we've parsed a data file
             has_pre_{false}, ///< True if the data contains pre-piracy, non-treatment rows
             has_piracy_{false}, ///< True if the data contains LR piracy treatment rows
             has_piracy_sr_{false}, ///< True if the data contains SR piracy treatment rows
             has_policy_{false}, ///< True if the data contains LR policy treatment rows
             has_policy_sr_{false}, ///< True if the data contains SR policy treatment rows
             require_pre_{false}, ///< True if pre data is required
             require_piracy_{false}, ///< True if piracy data is required
             require_policy_{false}, ///< True if policy data is required
             require_sr_{false}; ///< True if short-run data is required for each type of long-run data
        /// The number of treatment row observations per source data rows (i.e. per simulation)
        unsigned int rows_per_sim_{0};

        /// Field name to column map.
        std::unordered_map<std::string, unsigned> data_column_;

        /// Cache of weak pointers to returned SimpleVariables.
        mutable std::unordered_map<std::string, std::shared_ptr<SimpleVariable>> var_cache_;

        /** Generates a row of data from the CSVParser at its current position.
         *
         * \param csv the CSVParser positioned at the desired row.
         * \param newrow the matrix row in which values will be set
         * \param prefix the prefix (such as "pre." or "policy.SR.") to add to non-param. data
         * columns.
         */
        void generateRow(const CSVParser &csv,
                Eigen::Ref<Eigen::RowVectorXd, 0, Eigen::InnerStride<>> newrow,
                const std::string &prefix) const;

};

}}
