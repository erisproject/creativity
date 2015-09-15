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

struct tabulation_options {
    TableFormat format;
    unsigned precision;
    /** Constructs a tabulation_options struct.
     *
     * \param format the table format, one of the TableFormat enum values; defaults to Text
     * \param precision the precision of double values; defaults to 6
     */
    tabulation_options(TableFormat format = TableFormat::Text, unsigned precision = 6) : format{format}, precision{precision} {}
};


/** Function that takes an Eigen matrix and optional row and column names and displays the results.
 *
 * \param matrix the matrix (or matrix-like Eigen object) to display
 * \param options the tabulation_options controlling display output
 * \param rownames the row names to use.  If too short (including empty, the default), unspecified
 * rows will be named `[r,]` where r is the row index beginning at 0.
 * \param colnames the column names to use.  If too short (including empty, the default), unprovided
 * columns will be named `[,c]`, where c is the column index beginning at 0.
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
 * \param extracol is an extra column of string values to display to the right of the last matrix
 * column.  If empty (the default), no extra column is added.  The first value is the column title
 * (if any), the remaining `matrix.rows()` values are values to display.  If the matrix is too
 * short, blank entries will be used for missing values.
 */
std::string tabulate(
        const Equation &equation,
        const tabulation_options &options = tabulation_options{},
        const std::vector<std::string> &rownames = {},
        const std::vector<std::string> &extracol = {});

/** Function that takes an SUR object (i.e. of multiple equations) and prints them out, side-by-side
 * in a table; each y variable heading is prefixed with "Eq. n: ", with n starting from 1.
 *
 * \param sur the SUR model to display
 * \param options the tabulation_options controlling display output
 * \param rownames the row names to use.  If too short (including empty, the default), unspecified
 * rows will be named `[r,]` where r is the row index beginning at 0.
 * \param extracols is an vector of vectors of extra column of string values to display to the right of each model results in
 * result output; the first `sur.equations().size()` vectors will be displayed as columns, using the first `surthe last matrix
 * column.  If empty (the default), no extra column is added.  The first value is the column title
 * (if any), the remaining `matrix.rows()` values are values to display.  If the matrix is too
 * short, blank entries will be used for missing values.
 */
std::string tabulate(
        const SUR &sur,
        const tabulation_options &options = tabulation_options{},
        const std::vector<std::string> &rownames = {});


}}
