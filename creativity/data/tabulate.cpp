#include "creativity/data/tabulate.hpp"
#include <sstream>
#include <iomanip>

#include <iostream> // DEBUG

using namespace Eigen;

namespace creativity { namespace data {

// Extract the rowname, using '[r,]' if r exceeds the rownames length
inline std::string rowname(const std::vector<std::string> &rownames, unsigned r) {
    return (r >= rownames.size()) ? "[" + std::to_string(r) + ",]" : rownames[r];
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

    std::string interrow("  ");

    std::vector<unsigned> col_width(matrix.cols());
    for (unsigned c = 0; c < matrix.cols(); c++) {
        auto &width = col_width[c];
        width = colname(colnames, c).length();
        for (unsigned r = 0; r < matrix.rows(); r++) {
            std::ostringstream out;
            out.precision(options.precision);
            out << matrix(r,c);
            unsigned len = out.str().length();
            if (len > width) width = len;
        }
    }
    size_t extraw = 0;
    for (auto &s : extracol) extraw = std::max(extraw, s.length());

    unsigned rownames_width = 0;
    for (unsigned r = 0; r < matrix.rows(); r++) {
        unsigned len = rowname(rownames, r).length();
        if (len > rownames_width) rownames_width = len;
    }

    std::ostringstream table;
    table.precision(options.precision);
    table << options.indent << std::setw(rownames_width) << "";
    for (unsigned c = 0; c < matrix.cols(); c++) {
        table << interrow << std::setw(col_width[c]) << colname(colnames, c);
    }
    if (extraw > 0) table << std::left << std::setw(extraw) << extracol[0];
    table << "\n";

    for (unsigned r = 0; r < matrix.rows(); r++) {
        table << options.indent << std::setw(rownames_width) << rowname(rownames, r);
        for (unsigned c = 0; c < matrix.cols(); c++) {
            table << interrow << std::right << std::setw(col_width[c]) << matrix(r, c);
        }
        if (r+1 < extracol.size()) table << std::left << std::setw(extraw) << extracol[r+1];
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
}

std::string tabulate_latex(
        const Ref<const MatrixXd> &matrix,
        const tabulation_options &options,
        const std::vector<std::string> &rownames,
        const std::vector<std::string> &colnames,
        const std::vector<std::string> &extracol) {
    return "LaTeX FIXME\n";
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
    Eigen::MatrixXd mat(equation.depVar().size(), equation.numVars() + 1);
    std::vector<std::string> varnames;
    varnames.emplace_back(equation.depVar().name());
    equation.depVar().populate(mat.col(0));
    unsigned i = 1;
    for (const Variable &var : equation) {
        var.populate(mat.col(i++));
        varnames.emplace_back(var.name());
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
        std::cerr << eq.depVar().name() << ": " << eq.depVar().size() << " rows\n";
    }

    Eigen::MatrixXd mat(sur.equations().front().depVar().size(), cols);
    std::vector<std::string> varnames;
    unsigned eqnum = 1, i = 0;
    for (auto &eq : sur.equations()) {
        varnames.push_back("Eq. " + std::to_string(eqnum++) + ": " + eq.depVar().name());
        eq.depVar().populate(mat.col(i++));
        for (const Variable &var : eq) {
            varnames.emplace_back(var.name());
            var.populate(mat.col(i++));
        }
    }

    return tabulate(mat, options, rownames, varnames);
}

}}
