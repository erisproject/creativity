#include "creativity/data/CSVParser.hpp"
#include "creativity/data/SUR.hpp"
#include "creativity/data/tabulate.hpp"
#include "creativity/data/simdata.hpp"
#include "creativity/data/util.hpp"
#include "creativity/data/Treatment.hpp"
#include "creativity/data/TreatmentFilter.hpp"
#include "creativity/cmdargs/Results.hpp"
#include "creativity/Policy.hpp"
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
using namespace creativity;
using namespace creativity::data;
using namespace eris;
using namespace Eigen;
namespace fs = boost::filesystem;

int main(int argc, char *argv[]) {
    Results args;
    try {
        args.parse(argc, argv);
    } catch (const std::exception &e) {
        std::cerr << "Invalid argument(s): " << e.what() << std::endl;
        exit(2);
    }

    // rawdata holds all the data read from the file.
    Treatment rawdata;
    rawdata.requirePre();
    if (args.analysis.shortrun) rawdata.requireSR();
    try {
        rawdata.readCSV(CSVParser(args.input));
    }
    catch (const std::exception &e) {
        std::cerr << "Error: failed to read input file `" << args.input << "': " << e.what() << "\n\n";
        exit(1);
    }

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

    // The minimum value of books_written_pc that needs to be satisified to count as writing occuring
    // during a period:
    double writing_threshold = 0.2;
    // If we're not doing shortrun, filter out any shortrun rows:
    auto shortrun_filter = args.analysis.shortrun
        ? [](bool,bool,bool,bool) { return true; }
        : [](bool, bool, bool, bool short_run) { return not short_run; };

    // If filtering by policy, data is the filtered version of rawdata; otherwise it is just the rawdata itself.
    std::unique_ptr<TreatmentFilter> data_policy;
    if (args.analysis.policy_filter) {
        data_policy.reset(new TreatmentFilter(rawdata, [&](const TreatmentFilter::Properties &p) {
                return args.analysis.policy == Policy((uint32_t) p.pre->value("param.policy"));
                },
                shortrun_filter));
    }
    Treatment &data = args.analysis.policy_filter ? *data_policy : rawdata;

// Macro that is true if at least writing_threshold books were written and the mean initial price is
// a number
#define WRITING_AND_MARKET(var) ((var)->value("books_written_pc") >= writing_threshold and not std::isnan((var)->value("book_p0")))
// Macro that is true whenever WRITING_AND_MARKET is true and either the shortrun variable does not
// exist, or else it also satisfied the above.
#define WRITING_AND_MARKET_SRLR(var) (WRITING_AND_MARKET(var) and (not args.analysis.shortrun or not (var##_SR) or WRITING_AND_MARKET(var##_SR)))

    // Observations with writing in every stage:
    TreatmentFilter data_writing_always(data, [&](const TreatmentFilter::Properties &p) {
            // Only allow simulations that have writing in every stage contained in the data
            return (not p.pre or WRITING_AND_MARKET(p.pre))
                and (not p.piracy or WRITING_AND_MARKET_SRLR(p.piracy))
                and (not p.policy or WRITING_AND_MARKET_SRLR(p.policy));
        },
        shortrun_filter
    );

    if (args.output.list_filtered_write) {
        for (const auto &source : data_writing_always.sourceFiles()) {
            out << source << "\n";
        }
        return 0; // Done, exit.
    }

    // Observations with no writing under piracy, but with writing in policy (and pre)
    TreatmentFilter data_no_piracy_writing(data, [&](const TreatmentFilter::Properties &p) {
            // Require that this data actually has pre, piracy, and policy:
            if (not p.pre or not p.piracy or not p.policy) return false;
            // writing in pre and policy:
            return WRITING_AND_MARKET(p.pre)
                and not WRITING_AND_MARKET_SRLR(p.piracy)
                and WRITING_AND_MARKET_SRLR(p.policy);
        },
        shortrun_filter
    );
    // Observations with pre writing but no policy writing:
    TreatmentFilter data_no_post_writing(data, [&](const TreatmentFilter::Properties &p) {
            // Require that this data actually has pre, piracy, and policy:
            if (not p.pre or not p.piracy or not p.policy) return false;
            // Require writing in pre"
            return WRITING_AND_MARKET(p.pre)
                and not WRITING_AND_MARKET_SRLR(p.piracy)
                and not WRITING_AND_MARKET_SRLR(p.policy);
        },
        shortrun_filter
    );
    // Observations with no pre writing:
    TreatmentFilter data_no_pre_writing(data, [&](const TreatmentFilter::Properties &p) {
            return p.pre and not WRITING_AND_MARKET(p.pre);
        },
        shortrun_filter
    );
    // Observations with pre and piracy, but not policy writing:
    TreatmentFilter data_no_pol_writing(data, [&](const TreatmentFilter::Properties &p) {
            // Require that this data actually has pre, piracy, and policy:
            if (not p.pre or not p.piracy or not p.policy) return false;
            // Require writing in pre"
            return WRITING_AND_MARKET(p.pre)
                and WRITING_AND_MARKET_SRLR(p.piracy)
                and not WRITING_AND_MARKET_SRLR(p.policy);
        },
        shortrun_filter
    );

    if (args.output.list_filtered_nowrite) {
        for (auto &list : {data_no_pre_writing, data_no_piracy_writing, data_no_pol_writing, data_no_post_writing}) {
            for (const auto &source : list.sourceFiles()) {
                out << source << "\n";
            }
        }
        return 0; // Done, exit
    }


    tabulation_options tabopts(args.format.type, args.format.precision, true);
    tabopts.indent = "    "; // Has an effect in text mode only

    if (not args.output.no_preamble) {
        out << tabulate_preamble(args.format.type);
    }

    if (args.analysis.summary) {
        if (!args.latex_variables.empty()) {
            out << "% Data summary stats (" << args.latex_variables << "):\n";
            auto lv = [&out, &args](const std::string &suffix, unsigned int value, const std::string &comment = "") {
                out << "\\newcommand{\\" << args.latex_variables << suffix << "}{" << value << "}";
                if (!comment.empty()) out << " % " << comment;
                out << "\n";
            };
            lv("Total", rawdata.simulations(), "total simulations in data set");
            lv("RowsPerSim", rawdata.rowsPerSimulation(), "analysis rows per simulation");
            lv("Policy", data.simulations(), "total simulations matching the requested policy (if applicable)");
            lv("WritingAlways", data_writing_always.simulations(), "# sims with writing in pre, piracy, and policy");
            lv("NoPreWriting", data_no_pre_writing.simulations(), "# sims with no pre-piracy writing");
            lv("NoPostWriting", data_no_post_writing.simulations(), "# sims with pre-piracy but no piracy or policy writing");
            lv("NoPiracyWriting", data_no_piracy_writing.simulations(), "# sims with pre-piracy and policy writing, but no piracy writing");
            lv("NoPolicyWriting", data_no_pol_writing.simulations(), "# sums with pre-piracy and piracy writing, but not no policy writing");
        }
        else {
            out << tabulate_escape(std::string("Data summary:\n") +
                    "    " + std::to_string(rawdata.simulations()) + " total simulations (with " + std::to_string(rawdata.rowsPerSimulation()) + " data rows per simulation)\n" +
                    (args.analysis.policy_filter ? "    " + std::to_string(data.simulations()) + " simulations match the requested policy\n" : std::string()) +
                    "    " + std::to_string(data_writing_always.simulations()) + " simulations with non-zero # books written during each stage\n" +
                    "    " + std::to_string(data_no_pre_writing.simulations()) + " simulations with zero books written during pre-piracy stage\n" +
                    "    " + std::to_string(data_no_piracy_writing.simulations()) + " simulations with zero books written under piracy, but writing resuming under the policy\n" +
                    "    " + std::to_string(data_no_post_writing.simulations()) + " simulations with zero books written during piracy and no recovery under the policy\n" +
                    "    " + std::to_string(data_no_pol_writing.simulations()) + " simulations with writing under piracy, but no writing under the policy\n",
                    args.format.type);
        }
    }

    if (args.analysis.write_or_not) {
        if (not args.no_headings)
            out << "Write-vs-nowrite analysis:\n==========================\n";

        auto param_opts = tabopts, cor_opts = tabopts;
        param_opts.showpoint = false;
        param_opts.dot_align = false;
        cor_opts.matrix.diagonal = false;
        std::vector<std::string> params({
                "readers", "density", "reader_step_mean", "reader_creation_scale_range", "creation_fixed",
                "creation_time", "cost_market", "cost_unit", "cost_piracy", "piracy_link_proportion"});
        if (args.analysis.initial) {
            for (const std::string &i : {"prob_write", "l_min", "l_range", "p_min", "p_range", "prob_keep", "keep_price"})
                params.push_back("initial." + i);
        }

        if (not args.analysis.policy_filter or args.analysis.policy.publicSharing()) {
            params.push_back("policy_public_sharing_tax");
        }
        if (not args.analysis.policy_filter or args.analysis.policy.publicVoting()) {
            for (const std::string &p : {"tax", "votes"})
                params.push_back("policy_public_voting_" + p);
        }
        if (not args.analysis.policy_filter or args.analysis.policy.catchPirates()) {
            // FIXME: there are more parameters (I'm just not varying them yet)
            for (const std::string &p : {"tax", "fine[1]"})
                params.push_back("policy_catch_" + p);
        }

        std::vector<std::string> pre_fields;
        if (args.analysis.pre) {
            for (auto &f :{
                "net_u", "book_p0", "book_sales", "book_profit", "book_quality", "books_written_pc"})
                pre_fields.push_back(f);
        }

        std::vector<std::string> params_abbrev;
        std::regex word_re("([a-zA-Z0-9])[a-zA-Z0-9]+");
        for (const auto &p : params) params_abbrev.push_back(std::regex_replace(p, word_re, "$1"));
        enum : unsigned { f_mean, f_se, f_min, f_5, f_25, f_median, f_75, f_95, f_max, /* last: captures size: */ num_fields };
        std::vector<std::string> colnames({"Mean", "s.e.", "Min", "5th %", "25th %", "Median", "75th %", "95th %", "Max"});

        if (args.format.latex) {
            param_opts.escape = false; // We put some latex in the title, so need non-escaping
            // Unfortunately that means we have to manually escape the column names:
            for (auto &col : colnames) col = tabulate_escape(col, param_opts);
        }

        for (auto d : {&data_no_pre_writing, &data_no_post_writing, &data_no_piracy_writing, &data_no_pol_writing, &data_writing_always}) {
            param_opts.title = "Parameter values for " + std::to_string(d->simulations()) + " simulations " +
                (d == &data_writing_always ? "with writing in all stages" :
                 d == &data_no_pre_writing ? "without pre-piracy writing" :
                 d == &data_no_piracy_writing ? "with no piracy writing, but recovery under the policy" :
                 d == &data_no_post_writing ? "without piracy or policy writing" :
                 "with piracy writing but not policy writing") + ":";

            if (args.format.latex) {
                param_opts.title = "\\underline{\\normalsize{" + param_opts.title + "}\\nopagebreak}\\nopagebreak";
            }

            std::vector<std::string> local_params, local_params_abbrev;
            for (size_t i = 0; i < params.size(); i++) {
                const auto &p = params[i];
                // For the no-pre writing case, prune out piracy and policy variables:
                if (d == &data_no_pre_writing and (
                            p.find("piracy") != p.npos or p.find("policy") != p.npos))
                    continue;

                local_params.push_back(p);
                local_params_abbrev.push_back(params_abbrev[i]);
            }

            // Look at conditional means of parameters in periods with no activity vs parameters in
            // periods with activity

            MatrixXd results(local_params.size() + pre_fields.size(), (unsigned) num_fields);
            MatrixXd X(d->data().rows() / d->rowsPerSimulation(), local_params.size() + pre_fields.size());
            int x_row = 0;
            for (const auto &p : local_params) {
                VectorXd rawvals = (*d)["param." + p]->values();
                VectorXd paramvals = VectorXd::Map(rawvals.data(), rawvals.size()/d->rowsPerSimulation(), InnerStride<Dynamic>(d->rowsPerSimulation()));
                std::sort(paramvals.data(), paramvals.data() + paramvals.size());
                results(x_row, f_min) = quantile(paramvals, 0);
                results(x_row, f_5) = quantile(paramvals, .05);
                results(x_row, f_25) = quantile(paramvals, .25);
                results(x_row, f_median) = quantile(paramvals, .5);
                results(x_row, f_75) = quantile(paramvals, .75);
                results(x_row, f_95) = quantile(paramvals, .95);
                results(x_row, f_max) = quantile(paramvals, 1);
                X.col(x_row++) = paramvals;
            }
            for (const auto &p : pre_fields) {
                VectorXd rawvals = (*d)[p]->values();
                // Start at 0, increment by nobs_per_sim: that should keep us in the pre-sim rows
                VectorXd prevals = VectorXd::Map(rawvals.data(), rawvals.size()/d->rowsPerSimulation(), InnerStride<Dynamic>(d->rowsPerSimulation()));
                std::sort(prevals.data(), prevals.data() + prevals.size());
                results(x_row, f_min) = quantile(prevals, 0);
                results(x_row, f_5) = quantile(prevals, .05);
                results(x_row, f_25) = quantile(prevals, .25);
                results(x_row, f_median) = quantile(prevals, .5);
                results(x_row, f_75) = quantile(prevals, .75);
                results(x_row, f_95) = quantile(prevals, .95);
                results(x_row, f_max) = quantile(prevals, 1);
                X.col(x_row++) = prevals;
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
            row_names.insert(row_names.end(), local_params.begin(), local_params.end());
            for (const auto &pre : pre_fields) row_names.push_back("pre." + pre);

            for (auto &row : row_names) row = tabulate_escape(row, param_opts);

            out << tabulate(results, param_opts, row_names, colnames) << "\n";

            if (args.analysis.write_or_not_corrcov) {
                // Convert lower triangle of covariance matrix into correlation values:
                for (int c = 0; c < corr.cols(); c++) {
                    for (int r = c+1; r < corr.rows(); r++)
                        corr(r,c) /= sqrt(corr(c,c) * corr(r,r));
                }


                out << "Correlations (below diagonal), variance (diagonal), and covariance (above diagonal):\n" << tabulate(
                        corr.topLeftCorner(local_params.size(), local_params.size()), cor_opts, local_params, local_params_abbrev) << "\n";
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
                "books_written_pc",
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

            if (data_writing_always.hasPolicy()) {
                if (data_writing_always.hasPolicySR())
                    eq % (data_writing_always["policy"]*data_writing_always["SR"])
                        + data_writing_always["policy"]*data_writing_always["LR"];
                else
                    eq % data_writing_always["policy"];
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
                title << "\nEquation " << j+1 << ": " << eq << ":";
                avg_opts.title = title.str();

                out << tabulate(avg_effects.summary(j), avg_opts, avg_effects.varNames(j), {"Coefficient", "std.err.", "t-stat", "p-value"}, avg_effects.pStars(j));
            }
        }
        else {
            // Condensed form: all results in one table
            Eigen::MatrixXd condensed(avg_effects.equations().size(), avg_effects.equations().at(0).numVars());
            std::vector<std::string> depvars, stars;
#define LATEX_NAME(v, textit) {v, "$\\mathit{" textit "}$"}
#define LATEX_NAME_SUB(v, textit, sub) {v, "$\\mathit{" textit "}_{" sub "}$"}
#define LATEX_NAME_QUANTILE(v, t) LATEX_NAME(v, t), LATEX_NAME_SUB(v"_5th", t, "(Q=.05)"), LATEX_NAME_SUB(v"_median", t, "(Q=.50)"), LATEX_NAME_SUB(v"_95th", t, "(Q=.95)")
            std::unordered_map<std::string, std::string> latex_name({
                    LATEX_NAME_QUANTILE("net_u", "net\\ utility"),
                    LATEX_NAME("books_written_pc", "books\\ written"),
                    LATEX_NAME_QUANTILE("book_quality", "book\\ quality"),
                    LATEX_NAME_SUB("book_p0", "book\\ p", "t=0"),
                    LATEX_NAME("book_revenue", "book\\ revenue"),
                    LATEX_NAME("book_profit", "book\\ profit"),
                    LATEX_NAME_QUANTILE("book_author_scale", "book\\ author\\ } {\\alpha"),
                    LATEX_NAME_QUANTILE("book_author_effort", "book\\ author\\ } {\\ell"),
                    LATEX_NAME("books_bought", "books\\ bought/reader"),
                    LATEX_NAME("books_pirated", "books\\ pirated/reader"),
                    LATEX_NAME("books_public", "books\\ public/reader"),
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
                // If all the p-value stars are the same, just put it once; otherwise separate each
                // with /
                if (std::all_of(st.cbegin() + 1, st.cend(), [&st](const std::string &s) { return s == st[0]; }))
                    stars.push_back(st[0]);
                else {
                    std::ostringstream s;
                    s << st[0];
                    for (size_t i = 1; i < st.size(); i++) s << '/' << st[i];
                    stars.push_back(s.str());
                }
            }

            if (not args.no_headings)
                avg_opts.title = "Average effects:";
            avg_opts.dot_align = false;
            avg_opts.showpoint = true;

            std::vector<std::string> columns;
            columns.push_back("Constant");
            if (args.format.latex) avg_opts.escape = false;

            std::string delta(args.format.latex ? "$\\Delta$ " : "d ");

            if (data_writing_always.hasPiracySR()) columns.push_back(delta + "Piracy (SR)");
            columns.push_back(delta + "Piracy");

            if (data_writing_always.hasPolicySR()) columns.push_back(delta + "Policy (SR)");
            columns.push_back(delta + "Policy");

            out << tabulate(condensed, avg_opts, depvars, columns, stars);
        }
    }


    if (args.analysis.marginal) {
        if (not args.no_headings)
            out << "\n\n\nMarginal effects:\n=================\n";

        SUR marg_effects;
        for (auto &y : {"net_u", "books_written_pc", "book_quality", "book_p0", "book_revenue", "book_profit"}) {
            Equation eq(data_writing_always[y]);
            eq % 1;
            if (data_writing_always.hasPiracy()) eq % data_writing_always["piracy"];
            if (data_writing_always.hasPolicy()) eq % data_writing_always["policy"];
            for (auto &x : {"param.density", "param.cost_market", "param.creation_time", "param.creation_fixed", "param.cost_unit"}) {
                eq % data_writing_always[x];
                if (data_writing_always.hasPiracy()) eq % (data_writing_always["piracy"] * data_writing_always[x]);
                if (data_writing_always.hasPolicy()) eq % (data_writing_always["policy"] * data_writing_always[x]);
            }
            // These don't get included in "pre":
            for (auto &x : {"param.cost_piracy", "param.piracy_link_proportion"}) {
                if (data_writing_always.hasPiracy()) eq % (data_writing_always["piracy"] * data_writing_always[x]);
                if (data_writing_always.hasPolicy()) eq % (data_writing_always["policy"] * data_writing_always[x]);
            }
            // These don't get included in "pre" or "piracy":
            for (auto &x : {"param.policy_public_sharing_tax", "param.policy_public_voting_tax", "param.policy_public_voting_votes",
                    "param.policy_catch_tax", "param.policy_catch_fine[1]"}) {
                if (data_writing_always.hasPolicy()) eq % (data_writing_always["policy"] * data_writing_always[x]);
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

