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

using creativity::cmdargs::Results;
using namespace creativity::data;
using namespace eris;
using namespace Eigen;

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
#define WRITING_AND_MARKET(var) ((var)->value("books_written") >= writing_threshold and not std::isnan((var)->value("book_p0")))
    // Observations with writing in every (LR) stage:
    TreatmentFilter data_writing_always(data, [&writing_threshold](const TreatmentFilter::Properties &p) {
            // Only allow simulations that have writing in every LR stage contained in the data
            return (not p.pre or WRITING_AND_MARKET(p.pre))
                and (not p.piracy or WRITING_AND_MARKET(p.piracy))
                and (not p.public_sharing or WRITING_AND_MARKET(p.public_sharing));
    });
    // Observations with no writing in LR piracy, but with writing in LR public (and pre)
    TreatmentFilter data_no_piracy_writing(data, [&writing_threshold](const TreatmentFilter::Properties &p) {
            // Require that this data actually has pre, piracy, and public:
            if (not p.pre or not p.piracy or not p.public_sharing) return false;
            // writing in pre and LR public:
            return WRITING_AND_MARKET(p.pre)
                and not WRITING_AND_MARKET(p.piracy)
                and WRITING_AND_MARKET(p.public_sharing);
    });
    // Observations with pre writing but no public (LR) writing:
    TreatmentFilter data_no_post_writing(data, [&writing_threshold](const TreatmentFilter::Properties &p) {
            // Require that this data actually has pre, piracy, and public:
            if (not p.pre or not p.piracy or not p.public_sharing) return false;
            // Require writing in pre"
            return WRITING_AND_MARKET(p.pre)
                and not WRITING_AND_MARKET(p.piracy)
                and not WRITING_AND_MARKET(p.public_sharing);
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
                and WRITING_AND_MARKET(p.piracy)
                and not WRITING_AND_MARKET(p.public_sharing);
    });

    std::cout << "Data summary:\n" <<
        "    " << data.simulations() << " total simulations (with " << data.rowsPerSimulation() << " data rows per simulation)\n" <<
        "    " << data_writing_always.simulations() << " simulations with non-zero # books written during each stage\n" <<
        "    " << data_no_pre_writing.simulations() << " simulations with zero books written during pre-piracy stage\n" <<
        "    " << data_no_piracy_writing.simulations() << " simulations with zero books written under piracy, but writing resuming under public sharing\n" <<
        "    " << data_no_post_writing.simulations() << " simulations with zero books written during piracy and no recovery under public sharing\n" <<
        "    " << data_no_pub_writing.simulations() << " simulations with writing under piracy, but no writing under public sharing\n";

    if (args.analysis.write_or_not) {
        tabulation_options tab_opts(args.format.type, args.format.precision, "    ");
        tabulation_options cor_opts(args.format.type, args.format.precision, "    ");
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
        std::vector<std::string> colnames({
                "Mean", "s.e.", "Min", "5th %", "25th %", "Median", "75th %", "95th %", "Max",
                "Parameter"});

        for (auto d : {&data_writing_always, &data_no_piracy_writing, &data_no_pre_writing, &data_no_post_writing, &data_no_pub_writing}) {
            std::cout << "Parameter values for " << d->simulations() << " simulations " <<
                (d == &data_writing_always ? "with writing in all stages" :
                 d == &data_no_pre_writing ? "without pre-piracy writing" :
                 d == &data_no_piracy_writing ? "with no piracy writing, but recovery under public sharing" :
                 d == &data_no_post_writing ? "without piracy or public sharing writing" :
                 "with piracy writing but not public sharing writing") << "\n";
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

            std::vector<std::string> row_names(params);
            for (const auto &pre : pre_fields) row_names.push_back("pre." + pre);

            std::cout << tabulate(results, tab_opts, row_names, colnames) << "\n";

            if (args.analysis.write_or_not_corrcov) {
                // Convert lower triangle of covariance matrix into correlation values:
                for (int c = 0; c < corr.cols(); c++) {
                    for (int r = c+1; r < corr.rows(); r++)
                        corr(r,c) /= sqrt(corr(c,c) * corr(r,r));
                }


                std::cout << "Correlations (below diagonal) and covariance (above diagonal):\n" << tabulate(
                        corr.topLeftCorner(params.size(), params.size()), cor_opts, params, params_abbrev) << "\n";
            }
        }
    }


    SUR avg_effects;
    for (auto &y : {
            "net_u", "net_u_5th", "net_u_median", "net_u_95th",
            "books_written",
            "book_quality", "book_quality_5th", "book_quality_median", "book_quality_95th",
            "book_p0", "book_revenue", "book_profit",
            "book_author_scale", "book_author_scale_5th", "book_author_scale_median", "book_author_scale_95th",
            "book_author_effort", "book_author_effort_5th", "book_author_effort_median", "book_author_effort_95th"
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
                std::cerr << "Error: found NaN for " << avg_effects.equations()[i].depVar()->name() << ", source file " << data_writing_always.sourceFile(r) << "\n";
            }
        }
    }

    avg_effects.solve();
    std::cout << "Average effects:\n================\n" << avg_effects;



    SUR marg_effects;
    for (auto &y : {"net_u", "books_written", "book_quality", "book_p0", "book_revenue", "book_profit"}) {
        Equation eq(data[y]);
        eq % 1;
        if (data_writing_always.hasPiracy()) eq % data["piracy"];
        if (data_writing_always.hasPublic()) eq % data["public"];
        for (auto &x : {"param.density", "param.cost_market", "param.cost_unit", "param.creation_time", "param.creation_fixed"}) {
            eq % data[x];
            if (data_writing_always.hasPiracy()) eq % (data["piracy"] * data[x]);
            if (data_writing_always.hasPublic()) eq % (data["public"] * data[x]);
        }
        // These don't get included in "pre":
        for (auto &x : {"param.cost_piracy", "param.piracy_link_proportion"}) {
            if (data_writing_always.hasPiracy()) eq % (data["piracy"] * data[x]);
            if (data_writing_always.hasPublic()) eq % (data["public"] * data[x]);
        }
        // These don't get included in "pre" or "piracy":
        for (auto &x : {"param.public_sharing_tax"}) {
            if (data_writing_always.hasPublic()) eq % (data["public"] * data[x]);
        }
        marg_effects.add(eq);
    }

    marg_effects.solve();
    std::cout << "\n\n\nMarginal effects:\n=================\n" << marg_effects;
}

