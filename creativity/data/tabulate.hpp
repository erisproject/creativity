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
    TableFormat format = TableFormat::Text;
    /// The precision for numerical values
    unsigned precision = 6;
    /// Whether the decimal point and significant 0's are shown (see std::showpoint)
    bool showpoint = false;
    /// The indent (only applies when `format` is Text)
    std::string indent;
    /** The title of the table; if non-empty, this title will be displayed as a heading in some
     * format-specific way.  If empty (the default) no heading precedes the data.
     */
    std::string title;
    /** Whether to escape special characters (default) or leave them in formatted output.  This has
     * no effect for Text format, but matters for LaTeX (to escape things like %, _, etc.) and HTML
     * (<, &, etc.).
     */
    bool escape = true;
    /** If true, give row values and header left-alignment; false (the default) gives them right
     * alignment.
     */
    bool rows_align_left = false;
    /** If true (the default), give extracol values (if any) left-alignment; if false, use right
     * alignment.
     */
    bool extracol_align_left = true;
    /** For LaTeX and HTML, include a basic document preamble; if false (the default), output will
     * be a document fragment instead of a full document.
     *
     * For Text format, this has no effect.
     *
     * This is intended to allow for multiple tables to be displayed in the same (shell) document by
     * setting `preamble = true` for the first table, and `postamble = true` for the last.
     */
    bool preamble = false;
    /** For LaTeX and HTML, include a basic document postamble matching the preamble; if false (the
     * default), output will be a document fragment instead of a full document.
     *
     * For Text format, this has no effect.
     *
     * This is intended to allow for multiple tables to be displayed in the same (shell) document by
     * setting `preamble = true` for the first table, and `postamble = true` for the last.
     */
    bool postamble = false;
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
     * \param title the title of the table; defaults to empty.  Newline characters in the title will
     * be converted (if necessary) into newlines in the output format (e.g. `<br />` in HTML, `\\`
     * in LaTeX); other characters will not be escaped.
     * \param format the table format, one of the TableFormat enum values; defaults to Text
     * \param precision the precision of double values; defaults to 6
     * \param showpoint whether trailing 0s (and possibly the decimal point) should be shown for
     * displayed values; defaults to false (trailing 0s and/or decimal point are suppressed).
     */
    explicit tabulation_options(const std::string &title = "", TableFormat format = TableFormat::Text, unsigned precision = 6, bool showpoint = false)
        : format{format}, precision{precision}, showpoint{showpoint}, title{title} {}

    /** Constructs a tabulation_options struct without a title but with a specified format.  This is
     * exactly equivalent to calling the constructor with an empty string title plus the given
     * arguments.  Unlike the above, this constructor is not explicit, which allows implicit
     * conversion from a TableFormat to a tabulation_options struct.
     */
    tabulation_options(TableFormat format, unsigned precision = 6, bool showpoint = false)
        : tabulation_options("", format, precision, showpoint) {}
};


/** Function that takes an Eigen matrix and optional row and column names and displays the results.
 *
 * \param matrix the matrix (or matrix-like Eigen object) to display
 *
 * \param options the tabulation_options controlling display output
 *
 * \param rownames the row names to use.  The first value is the header for the rows; subsequent
 * values are the row names.  If this vector exactly equals the number of rows of the matrix, an
 * empty row name header (and so this is equivalent to passing a vector with an empty string
 * prepended to the beginning).  If the vector is smaller than the number of rows, the first value
 * is used as the row header, and unspecified rows will be named `[r,]` where r is the row index
 * beginning at 0.  If entirely empty, the row header will an the empty string.
 *
 * \param colnames the column names to use.  If too short (including empty, the default), unprovided
 * columns will be named `[,c]`, where c is the column index beginning at 0.
 *
 * \param extracol is an extra column of string values to display to the right of the last matrix
 * column.  If empty (the default), no extra column is added.  If the length of the vector exactly
 * equals the number of rows of the matrix, the extra column receives no header name; otherwise the
 * first value of the vector is used as the extra column header name.  If the vector is too short,
 * blank values will be used for any missing row values.  Values will be left-aligned in the column.
 */
std::string tabulate(
        const Eigen::Ref<const Eigen::MatrixXd> &matrix,
        const tabulation_options &options = tabulation_options{},
        const std::vector<std::string> &rownames = {},
        const std::vector<std::string> &colnames = {},
        const std::vector<std::string> &extracol = {});

/** Returns the preamble (if any) for the format contain within the given tabulation_options.  This
 * is called internally if tabulate() is called with the `preamble` option set to true, but can also
 * be called explicitly.
 *
 * \sa tabulation_options.preamble
 */
std::string tabulate_preamble(const tabulation_options &options);

/** Returns the postamble (if any) for the format contain within the given tabulation_options.  This
 * is called internally if tabulate() is called with the `postamble` option set to true, but can
 * also be called explicitly.
 *
 * \sa tabulation_options.postamble
 */
std::string tabulate_postamble(const tabulation_options &options);

/** Escapes a value according to the given tabulation options format.  This is primarily intended
 * for use when passing raw column or row values with the `escape` option set to false.
 *
 * The `options.escape` value is ignored for this method (that is, the method escapes even if the
 * `escape` value is set to false).
 *
 * \param in the input string to be escaped
 * \param options the tabulation options containing the `format` value specifying which type of
 * escaping is to be performed.  Can also be passed as simple a TableFormat value (using the
 * implicit tabulation_options constructor).
 * \returns a string with characters that need to be escaped replaced with suitable replacements for
 * the given `options.format` value.
 */
std::string tabulate_escape(const std::string &in, const tabulation_options &options);

/** Function that takes an Equation object and optional row names and displays the values of the
 * variables of the equation.  Column names come from the equation.
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

/** Function that takes an SUR object (i.e. of multiple equations) and returns a string of its data,
 * side-by-side in a table.  Each y variable heading is prefixed with "Eq. n: ", with n starting
 * from 1, then followed by the RHS variables of the equation, then the next y variable, etc.
 *
 * \param sur the SUR data to display
 * \param options the tabulation_options controlling display output
 * \param rownames the row names to use.  If too short (including empty, the default), unspecified
 * rows will be named `[r,]` where r is the row index beginning at 0.
 */
std::string tabulate(
        const SUR &sur,
        const tabulation_options &options = tabulation_options{},
        const std::vector<std::string> &rownames = {});

}}
