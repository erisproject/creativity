#include "creativity/data/tabulate.hpp"
#include <sstream>
#include <iomanip>
#include <regex>
#include <unordered_map>

using namespace Eigen;

namespace creativity { namespace data {

// Extract the rowname, using '[r,]' if r exceeds the rownames length.  If offset is not 0, then
// rownames[offset+i] corresponds to row i.
inline std::string rowname(const std::vector<std::string> &rownames, unsigned r, unsigned offset) {
    return (r+offset >= rownames.size())
        ? "[" + std::to_string(r) + ",]"
        : rownames[r+offset];
}
// Extract the colname, using '[,c]' if c exceeds the colnames length
inline std::string colname(const std::vector<std::string> &colnames, unsigned c) {
    return (c >= colnames.size()) ? "[," + std::to_string(c) + "]" : colnames[c];
}

std::string tabulate_text(
        const Ref<const MatrixXd> &matrix,
        const tabulation_options &options,
        const std::vector<std::string> &rownames,
        const std::vector<std::string> &colnames,
        const std::vector<std::string> &extracol) {

    std::string intercol("  ");

    std::vector<unsigned> col_width(matrix.cols());
    for (unsigned c = 0; c < matrix.cols(); c++) {
        auto &width = col_width[c];
        width = colname(colnames, c).length();
        for (unsigned r = 0; r < matrix.rows(); r++) {
            if (r < c ? options.matrix.upper : r > c ? options.matrix.lower : options.matrix.diagonal) {
                std::ostringstream out;
                out.precision(options.precision);
                out << matrix(r,c);
                unsigned len = out.str().length();
                if (len > width) width = len;
            }
        }
    }

    bool have_rowname_header = rownames.size() > 0 and rownames.size() != (unsigned) matrix.rows();
    unsigned rownames_width = have_rowname_header ? rownames[0].length() : 0;
    for (unsigned r = 0; r < matrix.rows(); r++) {
        unsigned len = rowname(rownames, r, have_rowname_header ? 1 : 0).length();
        if (len > rownames_width) rownames_width = len;
    }

    std::ostringstream table;
    // Set value precision and showpoint:
    table << std::setprecision(options.precision) << (options.showpoint ? std::showpoint : std::noshowpoint);

    // Output title, if there is one:
    if (not options.title.empty())
        table << options.title << "\n" << std::setfill('-') << std::setw(options.title.length()) << "" << std::setfill(' ') << "\n";

    // Output the header row:
    table << options.indent << std::setw(rownames_width) << (options.rows_align_left ? std::left : std::right) << (have_rowname_header ? rownames[0] : "");
    for (unsigned c = 0; c < matrix.cols(); c++)
        table << intercol << std::setw(col_width[c]) << colname(colnames, c);

    // Figure out if we need space for an extra column of strings:
    size_t extraw = 0;
    for (auto &s : extracol) extraw = std::max(extraw, s.length());
    unsigned extracol_offset = 0;
    if (extraw > 0) {
        std::string extracol_header;
        if (extracol.size() > 0 and extracol.size() != (unsigned) matrix.rows()) {
            extracol_offset = 1;
            extracol_header = extracol[0];
        }
        table << (options.extracol_align_left ? std::left : std::right) << std::setw(extraw) << extracol_header;
    }
    table << "\n";

    // Output matrix rows:
    for (unsigned r = 0; r < matrix.rows(); r++) {
        table << options.indent << std::setw(rownames_width) << (options.rows_align_left ? std::left : std::right)
            << rowname(rownames, r, have_rowname_header ? 1 : 0);
        for (unsigned c = 0; c < matrix.cols(); c++) {
            table << intercol << std::right << std::setw(col_width[c]);
            if (r < c ? options.matrix.upper : r > c ? options.matrix.lower : options.matrix.diagonal)
                table << matrix(r, c);
            else
                table << "";
        }
        // Append the extra column value, if applicable
        if (extraw > 0 and extracol.size() > r+extracol_offset)
            table << intercol << (options.extracol_align_left ? std::left : std::right) << std::setw(extraw) << extracol[r+extracol_offset];
        table << "\n";
    }

    return table.str();
}

std::string tabulate_html(
        const Ref<const MatrixXd> &matrix,
        const tabulation_options &options,
        const std::vector<std::string> &rownames,
        const std::vector<std::string> &colnames,
        const std::vector<std::string> &extracol) {
    return "HTML FIXME\n";
    if (false and (matrix.cols() or options.precision or rownames.empty() or colnames.empty() or extracol.empty())) {} // Silence unused var warnings
}


const std::unordered_map<char, std::string> latex_escape_map({
        {'#', "\\#"},
        {'$', "\\$"},
        {'%', "\\%"},
        {'&', "\\&"},
        {'_', "\\_"},
        {'{', "\\{"},
        {'}', "\\}"},
        {'<', "{\\textless}"},
        {'>', "{\\textgreater}"},
        {'^', "{\\textasciicircum}"},
        {'~', "{\\textasciitilde}"},
        {'\\', "{\\textbackslash}"},
        {'|', "{\\textbar}"},
        {'\n', "\\\\\n"}});

const std::unordered_map<char, std::string> html_escape_map({
        {'<', "&lt;"},
        {'>', "&gt;"},
        {'&', "&amp;"},
        {'"', "&quot;"},
        {'\n', "<br />\n"}});

// Escapes the given string using the escape characters/replacement pairs in the given map
std::string escape(const std::string &in, const std::unordered_map<char, std::string> &escape_map) {
    std::string escaped;
    auto end = escape_map.end();
    for (char c : in) {
        auto found = escape_map.find(c);
        if (found != end) escaped.append(found->second);
        else escaped.push_back(c);
    }
    return escaped;
}

// Like above, but only escapes when the escape option is specified, and uses the appropriate escape
// map based on the format option.
inline std::string escape(const std::string &in, const tabulation_options &opts) {
    if (opts.escape and opts.format == TableFormat::LaTeX) return escape(in, latex_escape_map);
    else if (opts.escape and opts.format == TableFormat::HTML) return escape(in, html_escape_map);
    else return in;
}

// The version advertised in the header (which always escapes).
std::string tabulate_escape(const std::string &in, const tabulation_options &opts) {
    if (opts.format == TableFormat::LaTeX) return escape(in, latex_escape_map);
    else if (opts.format == TableFormat::HTML) return escape(in, html_escape_map);
    else return in;
}

std::string tabulate_latex(
        const Ref<const MatrixXd> &matrix,
        const tabulation_options &options,
        const std::vector<std::string> &rownames,
        const std::vector<std::string> &colnames,
        const std::vector<std::string> &extracol) {

    std::ostringstream latex;

    if (options.preamble) latex << tabulate_preamble(options);

#define E(s) escape(s, options)

    latex << "\\begin{center}\n\n";
    if (not options.title.empty()) latex << std::regex_replace(E(options.title), std::regex("(?=\n)"), "\\\\") << "\n\n\\vspace{1em}\n\n";

    latex << "\\begin{tabular}{" << (options.rows_align_left ? "l" : "r");
    for (int i = 0; i < matrix.cols(); i++) latex << "r@{.}l";
    if (not extracol.empty()) latex << (options.extracol_align_left ? "l" : "r");
    latex << "}\n";

    if (rownames.size() > 0 and rownames.size() != (unsigned) matrix.rows()) latex << E(rownames[0]);
    latex << " &\n";

    for (unsigned c = 0; c < matrix.cols(); c++) {
        if (c > 0) latex << " &\n";
        latex << std::setw(2*c+2) << "" << "\\multicolumn{2}{c}{" << E(colname(colnames, c)) << "}";
    }
    latex << " \\\\[1ex]\n";


    std::regex dotexp(R"/((?:\.(\d*))?[eE](?:\+|(-))0*(\d+)$)/");
    std::string dotexp_replace("&$1$$\\times 10^{$2$3}$$");
    int rowname_offset = (rownames.size() > 0 and rownames.size() != (unsigned) matrix.rows()) ? 1 : 0;
    int extracol_offset = (extracol.size() > 0 and extracol.size() != (unsigned) matrix.rows()) ? 1 : 0;
    for (unsigned r = 0; r < matrix.rows(); r++) {
        latex << E(rowname(rownames, r, rowname_offset));
        for (unsigned c = 0; c < matrix.cols(); c++) {
            latex << " &\n" << std::setw(2*c+2) << "";

            if (r < c ? options.matrix.upper : r > c ? options.matrix.lower : options.matrix.diagonal) {
                double val = matrix(r, c);
                std::string value;
                {
                    std::ostringstream valstream;
                    valstream.precision(options.precision);
                    valstream << matrix(r, c);
                    value = valstream.str();
                }
                // If we have a finite value, replace a "." in the value with an & (the . gets put
                // in by the latex header specification).
                if (std::isfinite(val)) {
                    // If an floating-point notation number, use the dotexp regex to replace the .
                    // with a &, and replace the 'e' with '$\times 10^{...}$' (plus some minor
                    // cleanup in the regexp).
                    std::smatch m;
                    if (std::regex_search(value, m, dotexp)) {
                        latex << m.prefix() << m.format(dotexp_replace);
                    }
                    else {
                        // Not floating-point notation: if it has a ., replace it with an &;
                        // otherwise stick the & at the end.
                        auto dotpos = value.find_first_of('.');
                        if (dotpos == value.npos)
                            latex << value << "&";
                        else
                            latex << value.substr(0, dotpos) << "&" << value.substr(dotpos+1);
                    }
                }
                // Otherwise, for infinity, use the infinty symbol; otherwise (i.e. for NaN) just
                // leave the string "nan" as-is In either case, we make the field take up 2 columns
                // (since it has no "." separator)
                else {
                    latex << "\\multicolumn{2}{c}{";
                    if (std::isinf(val)) {
                        latex << "$";
                        if (val < 0) latex << "-";
                        latex << "\\infty$";
                    }
                    else {
                        latex << value;
                    }
                    latex << "}";
                }
            }
            else {
                latex << "\\multicolumn{2}{c}{}\n";
            }
        }

        if (not extracol.empty()) {
            latex << " &\n" << std::setw(2*matrix.cols()+2) << "";
            if (extracol.size() <= r + extracol_offset) latex << "{}";
            else {
                latex << E(extracol[extracol_offset + r]);
            }
        }

        latex << " \\\\\n";
    }

    latex << "\\end{tabular}\n\n\\end{center}\n";

    if (options.postamble) latex << tabulate_postamble(options);

#undef E

    return latex.str();
}

std::string tabulate_preamble(const tabulation_options &options) {
    if (options.format == TableFormat::LaTeX) return R"TeX(\documentclass[11pt]{article}
\usepackage[utf8]{inputenc}
\usepackage[left=2.54cm,top=2.54cm,right=2.54cm,bottom=2.54cm,nohead]{geometry}
\usepackage{longtable}

\begin{document}

\thispagestyle{empty}
)TeX";
    if (options.format == TableFormat::HTML) return R"HTML(<html>
<body>
)HTML";
    return "";
}

std::string tabulate_postamble(const tabulation_options &options) {
    if (options.format == TableFormat::LaTeX) return R"TeX(
\end{document}
)TeX";
    if (options.format == TableFormat::HTML) return R"HTML(
</body>
</html>
)HTML";
    return "";
}



std::string tabulate(
        const Ref<const MatrixXd> &matrix,
        const tabulation_options &options,
        const std::vector<std::string> &rownames,
        const std::vector<std::string> &colnames,
        const std::vector<std::string> &extracol) {

    switch (options.format) {
        case TableFormat::Text:
            return tabulate_text(matrix, options, rownames, colnames, extracol);
        case TableFormat::HTML:
            return tabulate_html(matrix, options, rownames, colnames, extracol);
        case TableFormat::LaTeX:
            return tabulate_latex(matrix, options, rownames, colnames, extracol);
    }
    return "";
}



std::string tabulate(
        const Equation &equation,
        const tabulation_options &options,
        const std::vector<std::string> &rownames) {
    Eigen::MatrixXd mat(equation.depVar()->size(), equation.numVars() + 1);
    std::vector<std::string> varnames;
    varnames.emplace_back(equation.depVar()->name());
    equation.depVar()->populate(mat.col(0));
    unsigned i = 1;
    for (const auto &var : equation) {
        var->populate(mat.col(i++));
        varnames.emplace_back(var->name());
    }

    return tabulate(mat, options, rownames, varnames);
}



std::string tabulate(
        const SUR &sur,
        const tabulation_options &options,
        const std::vector<std::string> &rownames) {
    unsigned cols = 0;
    for (auto &eq : sur.equations()) {
        cols += eq.numVars() + 1;
    }

    Eigen::MatrixXd mat(sur.equations().front().depVar()->size(), cols);
    std::vector<std::string> varnames;
    unsigned eqnum = 1, i = 0;
    for (auto &eq : sur.equations()) {
        varnames.push_back("Eq. " + std::to_string(eqnum++) + ": " + eq.depVar()->name());
        eq.depVar()->populate(mat.col(i++));
        for (const auto &var : eq) {
            varnames.emplace_back(var->name());
            var->populate(mat.col(i++));
        }
    }

    return tabulate(mat, options, rownames, varnames);
}

}}
