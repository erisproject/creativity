#include "creativity/data/CSVParser.hpp"
#include "creativity/data/SUR.hpp"
#include "creativity/data/tabulate.hpp"
#include "creativity/data/Data.hpp"
#include "creativity/data/Treatment.hpp"
#include "creativity/data/TreatmentFilter.hpp"
#include "creativity/cmdargs/Results.hpp"
#include <Eigen/Core>
#include <cerrno>
#include <exception>
#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <list>
#include <unordered_map>
#include <boost/filesystem/operations.hpp>

using creativity::cmdargs::Results;
using namespace creativity::data;
using namespace eris;
using namespace Eigen;
namespace fs = boost::filesystem;

int main(int argc, char *argv[]) {
    Results args;
    args.parse(argc, argv);

    // data holds all the data read from the file.
    Treatment data;
    try {
        data.readCSV(CSVParser(args.input));
    }
    catch (const std::exception &e) {
        std::cerr << "Error: failed to read input file `" << args.input << "': " << e.what() << "\n\n";
        exit(1);
    }

    // The minimum value of books_written that needs to be satisified to count as writing occuring
    // during a period:
    double writing_threshold = 0.2;
// Macro that is true if at least writing_threshold books were written and the mean initial price is
// a number
#define WRITING_AND_MARKET(var) ((var)->value("books_written") >= writing_threshold and not std::isnan((var)->value("book_p0")))
// Macro that is true whenever WRITING_AND_MARKET is true and either the shortrun variable does not
// exist, or else it also satisfied the above.
#define WRITING_AND_MARKET_SRLR(var) (WRITING_AND_MARKET(var) and (not (var##_SR) or WRITING_AND_MARKET(var##_SR)))
    // Observations with writing in every (LR) stage:
    TreatmentFilter data_writing_always(data, [&writing_threshold](const TreatmentFilter::Properties &p) {
            // Only allow simulations that have writing in every LR stage contained in the data
            return (not p.pre or WRITING_AND_MARKET(p.pre))
                and (not p.piracy or WRITING_AND_MARKET_SRLR(p.piracy))
                and (not p.public_sharing or WRITING_AND_MARKET_SRLR(p.public_sharing));
    });
    // Observations with no writing under piracy, but with writing in public (and pre)
    TreatmentFilter data_no_piracy_writing(data, [&writing_threshold](const TreatmentFilter::Properties &p) {
            // Require that this data actually has pre, piracy, and public:
            if (not p.pre or not p.piracy or not p.public_sharing) return false;
            // writing in pre and LR public:
            return WRITING_AND_MARKET(p.pre)
                and not WRITING_AND_MARKET_SRLR(p.piracy)
                and WRITING_AND_MARKET_SRLR(p.public_sharing);
    });
    // Observations with pre writing but no public (LR) writing:
    TreatmentFilter data_no_post_writing(data, [&writing_threshold](const TreatmentFilter::Properties &p) {
            // Require that this data actually has pre, piracy, and public:
            if (not p.pre or not p.piracy or not p.public_sharing) return false;
            // Require writing in pre"
            return WRITING_AND_MARKET(p.pre)
                and not WRITING_AND_MARKET_SRLR(p.piracy)
                and not WRITING_AND_MARKET_SRLR(p.public_sharing);
    });
    // Observations with no pre writing:
    TreatmentFilter data_no_pre_writing(data, [&writing_threshold](const TreatmentFilter::Properties &p) {
            return p.pre and not WRITING_AND_MARKET(p.pre);
    });
    // Observations with pre and piracy, but not public writing:
    TreatmentFilter data_no_pub_writing(data, [&writing_threshold](const TreatmentFilter::Properties &p) {
            // Require that this data actually has pre, piracy, and public:
            if (not p.pre or not p.piracy or not p.public_sharing) return false;
            // Require writing in pre"
            return WRITING_AND_MARKET(p.pre)
                and WRITING_AND_MARKET_SRLR(p.piracy)
                and not WRITING_AND_MARKET_SRLR(p.public_sharing);
    });

    std::ofstream f;
    if (not args.output.filename.empty()) {
        if (not args.output.overwrite and fs::exists(args.output.filename)) {
            std::cerr << "Error: `" << args.output.filename << "' already exists; specify a different file or add `--overwrite' option to overwrite\n";
            exit(1);
        }

        try {
            f.exceptions(f.failbit | f.badbit);
            f.open(args.output.filename, std::ios_base::trunc);
        }
        catch (std::ios_base::failure &c) {
            // Catch, wrap message, and rethrow:
            throw std::ios_base::failure("Unable to open " + args.output.filename + ": " + strerror(errno));
        }
    }
    std::ostream &out = args.output.filename.empty() ? std::cout : f;


    tabulation_options tabopts(args.format.type, args.format.precision, true);
    tabopts.indent = "    "; // Has an effect in text mode only

    if (not args.output.no_preamble) {
        out << tabulate_preamble(args.format.type);
    }

    if (args.analysis.summary) {
        out << tabulate_escape(std::string("Data summary:\n") +
                "    " + std::to_string(data.simulations()) + " total simulations (with " + std::to_string(data.rowsPerSimulation()) + " data rows per simulation)\n" +
                "    " + std::to_string(data_writing_always.simulations()) + " simulations with non-zero # books written during each stage\n" +
                "    " + std::to_string(data_no_pre_writing.simulations()) + " simulations with zero books written during pre-piracy stage\n" +
                "    " + std::to_string(data_no_piracy_writing.simulations()) + " simulations with zero books written under piracy, but writing resuming under public sharing\n" +
                "    " + std::to_string(data_no_post_writing.simulations()) + " simulations with zero books written during piracy and no recovery under public sharing\n" +
                "    " + std::to_string(data_no_pub_writing.simulations()) + " simulations with writing under piracy, but no writing under public sharing\n",
                args.format.type);
    }

    if (args.analysis.write_or_not) {
        if (not args.no_headings)
            out << "Write-vs-nowrite analysis:\n==========================\n";

        auto param_opts = tabopts, cor_opts = tabopts;
        param_opts.showpoint = false;
        param_opts.dot_align = false;
        cor_opts.matrix.diagonal = false;
        const std::vector<std::string> params({
                "readers", "density", "reader_step_mean", "reader_creation_scale_range", "creation_fixed",
                "creation_time", "cost_market", "cost_unit", "cost_piracy", "initial.prob_write", "initial.l_min",
                "initial.l_range", "initial.p_min", "initial.p_range", "initial.prob_keep", "initial.keep_price",
                "piracy_link_proportion", "public_sharing_tax"});
        const std::vector<std::string> pre_fields({
                "net_u", "book_p0", "book_sales", "book_profit", "book_quality", "books_written"});
        std::vector<std::string> params_abbrev;
        std::regex word_re("([a-zA-Z0-9])[a-zA-Z0-9]+");
        for (const auto &p : params) params_abbrev.push_back(std::regex_replace(p, word_re, "$1"));
        enum : unsigned { f_mean, f_se, f_min, f_5, f_25, f_median, f_75, f_95, f_max, /* last: captures size: */ num_fields };
        std::vector<std::string> colnames({"Mean", "s.e.", "Min", "5th %", "25th %", "Median", "75th %", "95th %", "Max"});

        for (auto d : {&data_writing_always, &data_no_piracy_writing, &data_no_pre_writing, &data_no_post_writing, &data_no_pub_writing}) {
            param_opts.title = "Parameter values for " + std::to_string(d->simulations()) + " simulations " +
                (d == &data_writing_always ? "with writing in all stages" :
                 d == &data_no_pre_writing ? "without pre-piracy writing" :
                 d == &data_no_piracy_writing ? "with no piracy writing, but recovery under public sharing" :
                 d == &data_no_post_writing ? "without piracy or public sharing writing" :
                 "with piracy writing but not public sharing writing") + ":";
            // Look at conditional means of parameters in periods with no activity vs parameters in
            // periods with activity

            MatrixXd results(params.size() + pre_fields.size(), (unsigned) num_fields);
            MatrixXd X(d->data().rows() / d->rowsPerSimulation(), params.size() + pre_fields.size());
            int i = 0;
            for (const auto &p : params) {
                VectorXd rawvals = (*d)["param." + p]->values();
                VectorXd paramvals = VectorXd::Map(rawvals.data(), rawvals.size()/d->rowsPerSimulation(), InnerStride<Dynamic>(d->rowsPerSimulation()));
                std::sort(paramvals.data(), paramvals.data() + paramvals.size());
                results(i, f_min) = quantile(paramvals, 0);
                results(i, f_5) = quantile(paramvals, .05);
                results(i, f_25) = quantile(paramvals, .25);
                results(i, f_median) = quantile(paramvals, .5);
                results(i, f_75) = quantile(paramvals, .75);
                results(i, f_95) = quantile(paramvals, .95);
                results(i, f_max) = quantile(paramvals, 1);
                X.col(i++) = paramvals;
            }
            for (const auto &p : pre_fields) {
                VectorXd rawvals = (*d)[p]->values();
                // Start at 0, increment by nobs_per_sim: that should keep us in the pre-sim rows
                VectorXd prevals = VectorXd::Map(rawvals.data(), rawvals.size()/d->rowsPerSimulation(), InnerStride<Dynamic>(d->rowsPerSimulation()));
                std::sort(prevals.data(), prevals.data() + prevals.size());
                results(i, f_min) = quantile(prevals, 0);
                results(i, f_5) = quantile(prevals, .05);
                results(i, f_25) = quantile(prevals, .25);
                results(i, f_median) = quantile(prevals, .5);
                results(i, f_75) = quantile(prevals, .75);
                results(i, f_95) = quantile(prevals, .95);
                results(i, f_max) = quantile(prevals, 1);
                X.col(i++) = prevals;
            }
            RowVectorXd mu = X.colwise().mean();
            results.col(f_mean) = mu.transpose();
            // Demean X:
            X.rowwise() -= mu;
            // Calculate covariance matrix:
            MatrixXd corr = (X.transpose() * X) / (X.rows()-1);
            // Extract the diagonals for the se values:
            results.col(f_se) = corr.diagonal().cwiseSqrt();

            std::vector<std::string> row_names;
            row_names.push_back("Parameter");
            row_names.insert(row_names.end(), params.begin(), params.end());
            for (const auto &pre : pre_fields) row_names.push_back("pre." + pre);

            out << tabulate(results, param_opts, row_names, colnames) << "\n";

            if (args.analysis.write_or_not_corrcov) {
                // Convert lower triangle of covariance matrix into correlation values:
                for (int c = 0; c < corr.cols(); c++) {
                    for (int r = c+1; r < corr.rows(); r++)
                        corr(r,c) /= sqrt(corr(c,c) * corr(r,r));
                }


                out << "Correlations (below diagonal), variance (diagonal), and covariance (above diagonal):\n" << tabulate(
                        corr.topLeftCorner(params.size(), params.size()), cor_opts, params, params_abbrev) << "\n";
            }
        }
    }


    if (args.analysis.average) {
        if (not args.no_headings and not args.condensed)
            out << "\n\n\nAverage effects:\n================\n";

        SUR avg_effects;
#define VALUE_PLUS_DIST(v) v, v"_5th", v"_median", v"_95th"
        for (auto &y : {
                VALUE_PLUS_DIST("net_u"),
                "books_written",
                VALUE_PLUS_DIST("book_quality"),
                "book_p0",
                "book_revenue",
                "book_profit",
                VALUE_PLUS_DIST("book_author_scale"),
                VALUE_PLUS_DIST("book_author_effort")
                }) {
            Equation eq(data_writing_always[y]);
            eq % 1;
            if (data_writing_always.hasPiracy()) {
                if (data_writing_always.hasPiracySR())
                    eq % (data_writing_always["piracy"]*data_writing_always["SR"])
                        + data_writing_always["piracy"]*data_writing_always["LR"];
                else
                    eq % data_writing_always["piracy"];
            }
            else if (data_writing_always.hasPiracy())
                eq % data_writing_always["piracy"];

            if (data_writing_always.hasPublic()) {
                if (data_writing_always.hasPublicSR())
                    eq % (data_writing_always["public"]*data_writing_always["SR"])
                        + data_writing_always["public"]*data_writing_always["LR"];
                else
                    eq % data_writing_always["public"];
            }

            avg_effects.add(std::move(eq));
        }
        avg_effects.gather();
        for (unsigned i = 0; i < avg_effects.equations().size(); i++) {
            auto yi = avg_effects.y(i);
            for (unsigned r = 0; r < yi.size(); r++) {
                if (std::isnan(yi[r])) {
                    std::cerr << "Error: found NaN for " << avg_effects.equations()[i].depVar()->name() << "\n"", source file " << data_writing_always.sourceFile(r) << "\n";
                }
            }
        }

        avg_effects.solve();

        auto avg_opts = tabopts;
        if (not args.condensed) {
            // Full set of tables
            for (unsigned j = 0; j < avg_effects.equations().size(); j++) {
                auto &eq = avg_effects.equations()[j];
                std::ostringstream title;
                title << "Equation " << j+1 << ": " << eq << ":";
                avg_opts.title = title.str();

                out << tabulate(avg_effects.summary(j), avg_opts, avg_effects.varNames(j), {"Coefficient", "std.err.", "t-stat", "p-value"}, avg_effects.pStars(j));
            }
        }
        else {
            // Condensed form: all results in one table
            Eigen::MatrixXd condensed(avg_effects.equations().size(), 3);
            std::vector<std::string> depvars, stars;
#define LATEX_NAME(v, textit) {v, "$\\mathit{" textit "}$"}
#define LATEX_NAME_SUB(v, textit, sub) {v, "$\\mathit{" textit "}_{" sub "}$"}
#define LATEX_NAME_QUANTILE(v, t) LATEX_NAME_SUB(v, t), LATEX_NAME_SUB(v"_5th", t, "(Q=.05)"), LATEX_NAME_SUB(v"_median", t, "(Q=.50)"), LATEX_NAME_SUB(v"_95th", t, "(Q=.95)")
            std::unordered_map<std::string, std::string> latex_name({
                    LATEX_NAME_QUANTILE("net_u", "net\\ utility"),
                    LATEX_NAME("books_written", "books\\ written"),
                    LATEX_NAME_QUANTILE("book_quality", "book\\ quality"),
                    LATEX_NAME_SUB("book_p0", "book\\ p", "t=0"),
                    LATEX_NAME("book_revenue", "book\\ revenue"),
                    LATEX_NAME("book_profit", "book\\ profit"),
                    LATEX_NAME_QUANTILE("book_author_scale", "book\\ author\\ } {\\alpha"),
                    LATEX_NAME_QUANTILE("book_author_effort", "book\\ author\\ } {\\ell")
                    });

            for (unsigned j = 0; j < avg_effects.equations().size(); j++) {
                std::string y_name = avg_effects.equations()[j].depVar()->name();
                if (args.format.latex) {
                    if (latex_name.count(y_name) > 0) {
                        depvars.push_back(latex_name.at(y_name));
                    }
                    else {
                        depvars.push_back(tabulate_escape(y_name, TableFormat::LaTeX));
                        std::cerr << "Warning: no LaTeX name substitute found for " << y_name << "\n";
                    }
                }
                else
                    depvars.push_back(y_name);

                condensed.row(j) = avg_effects.beta(j).transpose();
                auto st = avg_effects.pStars(j);
                if (st[0] == st[1] and st[0] == st[2]) stars.push_back(st[0]);
                else stars.push_back(st[0] + "/" + st[1] + "/" + st[2]);
            }

            if (not args.no_headings)
                avg_opts.title = "Average effects:";
            avg_opts.dot_align = false;
            avg_opts.showpoint = true;

            std::vector<std::string> columns;
            columns.push_back("Constant");
            if (args.format.latex) {
                avg_opts.escape = false;
                columns.push_back("$\\Delta$ Piracy");
                columns.push_back("$\\Delta$ Public");
            }
            else {
                columns.push_back("Piracy");
                columns.push_back("Public");
            }
            out << tabulate(condensed, avg_opts, depvars, columns, stars);
        }
    }


    if (args.analysis.marginal) {
        if (not args.no_headings)
            out << "\n\n\nMarginal effects:\n=================\n";

        SUR marg_effects;
        for (auto &y : {"net_u", "books_written", "book_quality", "book_p0", "book_revenue", "book_profit"}) {
            Equation eq(data_writing_always[y]);
            eq % 1;
            if (data_writing_always.hasPiracy()) eq % data_writing_always["piracy"];
            if (data_writing_always.hasPublic()) eq % data_writing_always["public"];
            for (auto &x : {"param.density", "param.cost_market", "param.creation_time", "param.creation_fixed", "param.cost_unit"}) {
                eq % data_writing_always[x];
                if (data_writing_always.hasPiracy()) eq % (data_writing_always["piracy"] * data_writing_always[x]);
                if (data_writing_always.hasPublic()) eq % (data_writing_always["public"] * data_writing_always[x]);
            }
            // These don't get included in "pre":
            for (auto &x : {"param.cost_piracy", "param.piracy_link_proportion"}) {
                if (data_writing_always.hasPiracy()) eq % (data_writing_always["piracy"] * data_writing_always[x]);
                if (data_writing_always.hasPublic()) eq % (data_writing_always["public"] * data_writing_always[x]);
            }
            // These don't get included in "pre" or "piracy":
            for (auto &x : {"param.public_sharing_tax"}) {
                if (data_writing_always.hasPublic()) eq % (data_writing_always["public"] * data_writing_always[x]);
            }
            marg_effects.add(eq);
        }

        marg_effects.solve();

        auto marg_opts = tabopts;
        for (unsigned j = 0; j < marg_effects.equations().size(); j++) {
            auto &eq = marg_effects.equations()[j];
            std::ostringstream title;
            title << "Equation " << j+1 << ": " << eq << ":";
            marg_opts.title = title.str();

            out << tabulate(marg_effects.summary(j), marg_opts, marg_effects.varNames(j), {"Coefficient", "std.err.", "t-stat", "p-value"}, marg_effects.pStars(j));
        }
    }

    if (not args.output.no_preamble) {
        out << tabulate_postamble(args.format.type);
    }
}

