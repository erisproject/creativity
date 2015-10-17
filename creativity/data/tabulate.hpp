#pragma once
#include "creativity/data/Equation.hpp"
#include "creativity/data/SUR.hpp"
#include <Eigen/Core>
#include <string>

namespace creativity { namespace data {

/// Enum for supported table output formats.
enum class TableFormat {
    Text, ///< Plain-text format
    HTML, ///< HTML table format
    LaTeX ///< LaTeX output
};

/// Struct holding various options controlling tabulation output
struct tabulation_options {
    /// The format (Text, HTML, or LaTeX)
    TableFormat format;
    /// The precision for numerical values
    unsigned precision;
    /// The indent (only applies when `format` is Text)
    std::string indent;
    /// The possible matrix sections to display
    struct {
        /// If true (the default), display the diagonal
        bool diagonal = true;
        /// If true (the default), display the part below the diagonal
        bool lower = true;
        /// If true (the default), display the part above the diagonal
        bool upper = true;
    } matrix;
    /** Constructs a tabulation_options struct.
     *
     * \param format the table format, one of the TableFormat enum values; defaults to Text
     * \param precision the precision of double values; defaults to 6
     * \param indent the per-line indent to use (only for `format == Text`); defaults to an empty
     * string
     */
    tabulation_options(TableFormat format = TableFormat::Text, unsigned precision = 6, const std::string &indent = "")
        : format{format}, precision{precision}, indent{indent} {}
};


/** Function that takes an Eigen matrix and optional row and column names and displays the results.
 *
 * \param matrix the matrix (or matrix-like Eigen object) to display
 * \param options the tabulation_options controlling display output
 * \param rownames the row names to use.  If too short (including empty, the default), unspecified
 * rows will be named `[r,]` where r is the row index beginning at 0.
 * \param colnames the column names to use.  If too short (including empty, the default), unprovided
 * columns will be named `[,c]`, where c is the column index beginning at 0.  If the size is exactly
 * equal to the matrix columns plus 1, the last value is used as a header for the row names.
 * \param extracol is an extra column of string values to display to the right of the last matrix
 * column.  If empty (the default), no extra column is added.  The first value is the column title
 * (if any), the remaining `matrix.rows()` values are values to display.  If the matrix is too
 * short, blank entries will be used for missing values.
 */
std::string tabulate(
        const Eigen::Ref<const Eigen::MatrixXd> &matrix,
        const tabulation_options &options = tabulation_options{},
        const std::vector<std::string> &rownames = {},
        const std::vector<std::string> &colnames = {},
        const std::vector<std::string> &extracol = {});

/** Function that takes an Equation object and optional row names and displays the results.  Column
 * names come from the equation.
 *
 * \param equation the equation object of values to display
 * \param options the tabulation_options controlling display output
 * \param rownames the row names to use.  If too short (including empty, the default), unspecified
 * rows will be named `[r,]` where r is the row index beginning at 0.
 */
std::string tabulate(
        const Equation &equation,
        const tabulation_options &options = tabulation_options{},
        const std::vector<std::string> &rownames = {});

/** Function that takes an SUR object (i.e. of multiple equations) and returns a string of them,
 * side-by-side in a table; each y variable heading is prefixed with "Eq. n: ", with n starting from
 * 1.
 *
 * \param sur the SUR model to display
 * \param options the tabulation_options controlling display output
 * \param rownames the row names to use.  If too short (including empty, the default), unspecified
 * rows will be named `[r,]` where r is the row index beginning at 0.
 */
std::string tabulate(
        const SUR &sur,
        const tabulation_options &options = tabulation_options{},
        const std::vector<std::string> &rownames = {});


}}
