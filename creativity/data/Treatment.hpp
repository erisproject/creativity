#pragma once
#include <creativity/data/CSVParser.hpp>
#include <creativity/data/Variable.hpp>
#include <Eigen/Core>
#include <memory>

namespace creativity { namespace data {

/** This class converts raw simulation data into multiple data rows, with piracy/public, SR/LR as
 * treatment effects on the base (pre-piracy, pre-public) row.  It additionally supports filtering
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

        /// True if this treatment data pre-piracy observations
        const bool& hasPre() const { return has_pre_; }

        /// True if the source data has piracy data
        const bool& hasPiracy() const { return has_piracy_; }

        /// True if the source data has public data
        const bool& hasPublic() const { return has_public_; }

        /// True if the source data has short-run piracy data
        const bool& hasPiracySR() const { return has_piracy_sr_; }

        /// True if the source data has short-run public data
        const bool& hasPublicSR() const { return has_public_sr_; }

        /// The number of data rows per observation.  This equals 1 plus the number of treatments.
        const unsigned int& rowsPerObservation() const { return rows_per_obs_; }

        /** The full set of data.
         *
         * \sa TreatmentFilter for extracting conditional data.
         */
        const Eigen::MatrixXd& data() const { return data_; }

        /** The simulation source values.  Data in row `x` is from source file at index `x /
         * rowsPerObservation()`.
         */
        const std::vector<std::string>& sourceFiles() const { return source_; }

        /** Return the source file associated with row r.  This is simply an alias for
         * `sourceFiles().at(r / rowsPerObservation())`.
         */
        const std::string& sourceFile(unsigned r) const { return sourceFiles().at(r / rowsPerObservation()); }

        /** Accesses the map of field name to column indices.  Data fields (e.g. "books_written")
         * are not prefixed with "pre.", "piracy.", etc.; parameters (e.g. "param.readers") are
         * prefixed; there are also dummies "piracy", "public", "SR", and "LR".
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

        bool has_pre_{false}, ///< True if the data contains pre-piracy, non-treatment rows
             has_piracy_{false}, ///< True if the data contains LR piracy treatment rows
             has_piracy_sr_{false}, ///< True if the data contains SR piracy treatment rows
             has_public_{false}, ///< True if the data contains LR public treatment rows
             has_public_sr_{false}; ///< True if the data contains SR public treatment rows
        /// The number of treatment row observations per source data rows (i.e. per simulation)
        unsigned int rows_per_obs_{0};

        /// Field name to column map.
        std::unordered_map<std::string, unsigned> data_column_;

        /// Cache of weak pointers to returned SimpleVariables.
        mutable std::unordered_map<std::string, std::shared_ptr<SimpleVariable>> var_cache_;

        /** Generates a row of data from the CSVParser at its current position.
         *
         * \param csv the CSVParser positioned at the desired row.
         * \param newrow the matrix row in which values will be set
         * \param prefix the prefix (such as "pre." or "public.SR.") to add to non-param. data
         * columns.
         */
        void generateRow(const CSVParser &csv,
                Eigen::Ref<Eigen::RowVectorXd, 0, Eigen::InnerStride<>> newrow,
                const std::string &prefix) const;

};

}}
