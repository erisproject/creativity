#include "creativity/data/CSVParser.hpp"
#include "creativity/data/SUR.hpp"
#include "creativity/data/tabulate.hpp"
#include "creativity/data/Data.hpp"
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

// rawdata holds the data for rows with positive number of books written under piracy/public,
// raw_no_writing holds the data for rows with books in pre but no books in one or more later situations.
// raw_no_writing_pir holds data for rows with no piracy writing, but with writing under public sharing
// raw_no_writing_pre holds the data for rows with no books written under pre-piracy.
MatrixXd rawdata, raw_no_writing, raw_no_writing_pre, raw_no_writing_pir;
std::unordered_map<std::string, std::shared_ptr<const SimpleVariable>>
    data_, data_nw_, data_nw_pir_, data_nw_pre_;
// Map row numbers to filenames:
std::vector<std::string> data_source, data_nw_source, data_nw_pir_source, data_nw_pre_source;
// Track fields that don't exist in some of the stages:
std::unordered_map<std::string, int> pre_nan_, piracy_nan_, public_nan_;
bool piracy_data = false, // True if we have piracy data
     public_data = false, // True if we have public data
     piracy_sr = false, // True if we have short-run piracy data
     public_sr = false; // True if we have short-run public sharing data

// Stupid wrapper class so that data["foo"] works (so data.at("foo") isn't needed)
class DataWrapper {
    public:
        DataWrapper(decltype(data_) &data) : d_{data} {}
        std::shared_ptr<const SimpleVariable> operator[](const std::string &field) {
            try { return d_.at(field); }
            catch (const std::out_of_range &re) {
                std::cerr << "Data exception: data field '" << field << "' does not exist\n";
                throw;
            }
        }
        const std::shared_ptr<const SimpleVariable>& operator[](const std::string &field) const {
            try { return d_.at(field); }
            catch (const std::out_of_range &re) {
                std::cerr << "Data exception: data field '" << field << "' does not exist\n";
                throw;
            }
        }
        bool has(const std::string &field) const { return d_.count(field) > 0; }
    private:
        decltype(data_) &d_;
};
DataWrapper data(data_), data_nw(data_nw_), data_nw_pre(data_nw_pre_), data_nw_pir(data_nw_pir_);


// Generate a row of data from a CSV row
void generateRow(
        const CSVParser &csv, // the CSV parser (current row will be used)
        Ref<RowVectorXd, 0, InnerStride<>> newrow, // Row to set
        const std::unordered_map<std::string, int> &data_column, // field to column map
        const std::unordered_map<std::string, int> &nans, // Fields to set to nan
        bool &no_writing, // Will be set true if books_written is less than 0.2 (won't be changed otherwise)
        const std::string &prefix, // Prefix to match (in addition to param.)
        const std::string &notprefix = "" // If non-empty, prefix to avoid even if prefix matches--typically `whatever.SR.`
        ) {

    for (unsigned fc = 0; fc < csv.fields().size(); fc++) {
        auto &d = csv.fields()[fc];
        if (
                d.substr(0, 6) == "param." or // Always include parameters
                (d.substr(0, prefix.size()) == prefix // Include things with the requested prefix
                    // ... unless `notprefix` is given and it matches the field
                    and not (!notprefix.empty() and d.substr(0, notprefix.size()) == notprefix))
           ) {
            newrow[data_column.at(d)] = csv.row()[fc];
        }
    }
    for (auto &nan : nans) newrow[nan.second] = std::numeric_limits<double>::quiet_NaN();
    if (piracy_data) newrow[data_column.at("piracy")] = prefix.substr(0, 7) == "piracy." ? 1 : 0;
    if (public_data) newrow[data_column.at("public")] = prefix.substr(0, 7) == "public." ? 1 : 0;
    if (piracy_sr or public_sr) {
        bool longrun = prefix.find(".SR.") == std::string::npos;
        newrow[data_column.at("LR")] = longrun;
        newrow[data_column.at("SR")] = !longrun;
    }

    if (
            (data_column.count(".books_written") and newrow[data_column.at(".books_written")] < 0.2)
        or
            (data_column.count(".book_p0") and std::isnan(newrow[data_column.at(".book_p0")]))
       )
        no_writing = true;
}



void readCSV(const std::string &filename) {
    CSVParser csv(filename);

    // The data contains values like pre.whatever, piracy.whatever, piracy.SR.whatever,
    // public.whatever, public.SR.whatever.  We need to convert those into five rows:
    // - one with `whatever' set to pre.whatever and piracy=public=SR=0, LR=1
    // - one with `whatever' set to piracy.SR.whatever, piracy=1, public=0, SR=1, LR=0
    // - one with `whatever' set to piracy.whatever, piracy=1, public=0, SR=0, LR=1
    // - one with `whatever' set to public.SR.whatever, piracy=0, public=1, SR=1, LR=0
    // - one with `whatever' set to public.whatever, piracy=0, public=1, SR=0, LR=1
    //
    // The data might not, however, have any piracy and/or public and/or short run data, in which
    // case we omit the relevant row(s).

    // This maps every CSV field (by name) into a column.  The mapping isn't unique, however: any
    // pre.whatever, piracy.whatever, and public.whatever should all map to the same column, albeit
    // with different rows, as described above.
    //
    // For the various pre/piracy/public fields, we also store the common location in `_whatever`,
    // so that it's easier to look up.  Note: elements with leading single underscores should be
    // considered special!
    std::unordered_map<std::string, int> data_column;

    int next_col = 0;
    for (auto &f : csv.fields()) {
        if (f[0] == '.') throw std::runtime_error("CSV file has invalid field `" + f + "': fields may not start with .");
        std::string data_field;
        if (f.substr(0, 4) == "pre.") data_field = f.substr(3);
        else if (f.substr(0, 10) == "piracy.SR.") { data_field = f.substr(9); piracy_sr = true; }
        else if (f.substr(0, 10) == "public.SR.") { data_field = f.substr(9); public_sr = true; }
        else if (f.substr(0, 7) == "piracy.") { data_field = f.substr(6); piracy_data = true; }
        else if (f.substr(0, 7) == "public.") { data_field = f.substr(6); public_data = true; }

        if (not data_field.empty()) {
            if (data_column.count(data_field)) {
                // One of the pre/piracy/public equivalent values for this value was already seen, reuse that column
                data_column.insert({f, data_column[data_field]});
            }
            else {
                // This field hasn't been seen before
                data_column.insert({f, next_col});
                data_column.insert({data_field, next_col});
                next_col++;
            }
        }
        else {
            // Not a pre/piracy/public field, so just copy it as-is
            data_column.insert({f, next_col});
            next_col++;
        }
    }

    if (piracy_data) data_column.insert({"piracy", next_col++}); // piracy dummy
    if (public_data) data_column.insert({"public", next_col++}); // public sharing dummy
    if (piracy_sr or public_sr) {
        data_column.insert({"SR", next_col++}); // short-run stage dummy
        data_column.insert({"LR", next_col++}); // long-run stage dummy
    }

    // Increase the matrix by 2520 row increments: 2520 is the lowest integer divisible by all
    // numbers from 1 to 10, so this will work properly for adding anywhere from 1 to 10 rows of
    // data per input row.  (Currently we have 1-5 possible rows, and 2520 is a reasonably sized
    // number for this, as well).
    constexpr unsigned rowincr = 2*2*2*3*3*5*7;
    rawdata.resize(0, next_col);
    raw_no_writing.resize(0, next_col);
    raw_no_writing_pre.resize(0, next_col);
    raw_no_writing_pir.resize(0, next_col);

    // Some fields don't have all three (for example, books_pirated has piracy_ and public_ values
    // but not a pre_ value; books_public_copies is only under public_), so we need to fill in some
    // NaNs for the ones where they don't exist: track the columns needing NaNs here:
    for (auto &pair : data_column) {
        if (pair.first[0] == '.') {
            if (not data_column.count("pre" + pair.first))
                pre_nan_.emplace(pair.first.substr(1), pair.second);
            if (piracy_data and not data_column.count("piracy" + pair.first))
                piracy_nan_.emplace(pair.first.substr(1), pair.second);
            if (public_data and not data_column.count("public" + pair.first))
                public_nan_.emplace(pair.first.substr(1), pair.second);
        }
    }

    unsigned raw_rows = 0, rawnw_rows = 0, rawnwpre_rows = 0, rawnwpir_rows = 0;
    while (csv.readRow()) {
        MatrixXd newrows(1 + piracy_data + public_data + piracy_sr + public_sr, rawdata.cols());
        unsigned rownum = 0;
        // Will be set to true if we find one or more periods with no writing/private marketing
        bool no_writing = false;

        // First row: pre values
        generateRow(csv, newrows.row(rownum++), data_column, pre_nan_, no_writing, "pre.");
        bool no_pre_writing = no_writing;

        bool no_pir_writing = false;
        // Second: short run piracy
        if (piracy_sr) {
            generateRow(csv, newrows.row(rownum++), data_column, piracy_nan_, no_writing, "piracy.SR.");
        }
        // Third: long run piracy
        if (piracy_data) {
            generateRow(csv, newrows.row(rownum++), data_column, piracy_nan_, no_pir_writing, "piracy.", "piracy.SR.");
        }
        // Fourth: short run public
        if (public_sr) {
            generateRow(csv, newrows.row(rownum++), data_column, public_nan_, no_writing, "public.SR.");
        }
        // Fifth: long run public
        if (public_data) {
            generateRow(csv, newrows.row(rownum++), data_column, public_nan_, no_writing, "public.", "public.SR.");
        }

        MatrixXd &rawmatrix = (no_pre_writing ? raw_no_writing_pre : no_writing ? raw_no_writing : no_pir_writing ? raw_no_writing_pir : rawdata);
        auto &rawrows = (no_pre_writing ? rawnwpre_rows : no_writing ? rawnw_rows : no_pir_writing ? rawnwpir_rows : raw_rows);
        auto &source = (no_pre_writing ? data_nw_pre_source : no_writing ? data_nw_source : no_pir_writing ? data_nw_pir_source : data_source);
        if (rawrows + newrows.rows() >= rawmatrix.rows()) {
            rawmatrix.conservativeResize(rawmatrix.rows() + rowincr, NoChange);
        }
        rawmatrix.middleRows(rawrows, newrows.rows()) = newrows;
        unsigned need = rawrows + newrows.rows();
        if (source.size() < need) source.resize(need);
        for (unsigned i = 0; i < newrows.rows(); i++) source[rawrows + i] = csv.rowSkipped().at("source");
        rawrows += newrows.rows();
    }

    // Resize the probably-oversized data matrices back to the number of rows we actually filled
    if (rawdata.rows() > raw_rows) rawdata.conservativeResize(raw_rows, NoChange);
    if (raw_no_writing.rows() > rawnw_rows) raw_no_writing.conservativeResize(rawnw_rows, NoChange);
    if (raw_no_writing_pre.rows() > rawnwpre_rows) raw_no_writing_pre.conservativeResize(rawnwpre_rows, NoChange);
    if (raw_no_writing_pir.rows() > rawnwpir_rows) raw_no_writing_pir.conservativeResize(rawnwpir_rows, NoChange);

    for (auto &dc : data_column) {
        if (dc.first.substr(0, 4) == "pre." or dc.first.substr(0, 7) == "public." or dc.first.substr(0, 7) == "piracy.")
            continue;
        std::string name(dc.first);
        if (dc.first[0] == '.') name = dc.first.substr(1);
        data_.insert({name, SimpleVariable::create(name, rawdata.col(dc.second))});
        data_nw_.insert({name, SimpleVariable::create(name, raw_no_writing.col(dc.second))});
        data_nw_pre_.insert({name, SimpleVariable::create(name, raw_no_writing_pre.col(dc.second))});
        data_nw_pir_.insert({name, SimpleVariable::create(name, raw_no_writing_pir.col(dc.second))});
    }
}

int main(int argc, char *argv[]) {
    Results args;
    args.parse(argc, argv);

    try {
        readCSV(args.input);
    }
    catch (const std::exception &e) {
        std::cerr << "Error: failed to read input file `" << args.input << "': " << e.what() << "\n\n";
    }

    unsigned nobs_per_sim = 1 + piracy_data + public_data + public_sr + piracy_sr;
    unsigned nobs_w_writing = data["net_u"]->size(),
             nobs_wo_writing = data_nw["net_u"]->size(),
             nobs_wo_pir_writing = data_nw_pir["net_u"]->size(),
             nobs_wo_pre_writing = data_nw_pre["net_u"]->size();
    unsigned sims_w_writing = nobs_w_writing / nobs_per_sim,
             sims_wo_writing = nobs_wo_writing / nobs_per_sim,
             sims_wo_pir_writing = nobs_wo_pir_writing / nobs_per_sim,
             sims_wo_pre_writing = nobs_wo_pre_writing / nobs_per_sim;
    std::cout << "Data summary:\n" <<
        "    " << nobs_w_writing + nobs_wo_writing + nobs_wo_pre_writing + nobs_wo_pir_writing <<
                    " total observations (from " << sims_w_writing+sims_wo_writing+sims_wo_pre_writing+sims_wo_pir_writing << " simulations)\n" <<
        "    " << nobs_w_writing << " observations (from " << sims_w_writing << " simulations) with non-zero # books written during each stage\n" <<
        "    " << nobs_wo_pre_writing << " observations (from " << sims_wo_pre_writing << " simulations) with zero books written during pre-piracy stage\n" <<
        "    " << nobs_wo_pir_writing << " observations (from " << sims_wo_pir_writing << " simulations) with zero books written during piracy, but positive books written during public sharing\n" <<
        "    " << nobs_wo_writing << " observations (from " << sims_wo_writing << " simulations) with zero books written during one or more situations\n";

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
        for (auto d : {&data, &data_nw_pre, &data_nw, &data_nw_pir}) {
            std::cout << "Parameter values for ";
            if (d == &data) std::cout << sims_w_writing << " simulations with writing in all stages:\n";
            else if (d == &data_nw_pre) std::cout << sims_wo_pre_writing << " simulations WITHOUT writing in last pre-piracy stages:\n";
            else if (d == &data_nw_pir) std::cout << sims_wo_pir_writing << " simulations with no writing in piracy stage, but with writing in public sharing stage:\n";
            else std::cout << sims_wo_writing << " simulations WITHOUT writing in either piracy or public sharing stages:\n";
            // Look at conditional means of parameters in periods with no activity vs parameters in
            // periods with activity

            MatrixXd results(params.size() + pre_fields.size(), (unsigned) num_fields);
            MatrixXd X((*d)["param.readers"]->size() / nobs_per_sim, params.size() + pre_fields.size());
            int i = 0;
            for (const auto &p : params) {
                VectorXd rawvals = (*d)[std::string("param.") + p]->values();
                VectorXd paramvals = VectorXd::Map(rawvals.data(), rawvals.size()/nobs_per_sim, InnerStride<Dynamic>(nobs_per_sim));
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
                VectorXd prevals = VectorXd::Map(rawvals.data(), rawvals.size()/nobs_per_sim, InnerStride<Dynamic>(nobs_per_sim));
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
        Equation eq(data[y]);
        eq % 1;
        if (piracy_sr)        eq % (data["piracy"]*data["SR"]) + data["piracy"]*data["LR"];
        else if (piracy_data) eq % data["piracy"];

        if (public_sr)        eq % (data["public"]*data["SR"]) + data["public"]*data["LR"];
        else if (public_data) eq % data["public"];

        avg_effects.add(std::move(eq));
    }
    avg_effects.gather();
    for (unsigned i = 0; i < avg_effects.equations().size(); i++) {
        auto yi = avg_effects.y(i);
        for (unsigned r = 0; r < yi.size(); r++) {
            if (std::isnan(yi[r])) {
                std::cerr << "Error: found NaN for " << avg_effects.equations()[i].depVar()->name() << ", source file " << data_source[r] << "\n";
            }
        }
    }

    avg_effects.solve();
    std::cout << "Average effects:\n================\n" << avg_effects;



    SUR marg_effects;
    for (auto &y : {"net_u", "books_written", "book_quality", "book_p0", "book_revenue", "book_profit"}) {
        Equation eq(data[y]);
        eq % 1;
        if (piracy_data) eq % data["piracy"];
        if (public_data) eq % data["public"];
        for (auto &x : {"param.density", "param.cost_market", "param.cost_unit", "param.creation_time", "param.creation_fixed"}) {
            if (not data.has(x)) continue;
            if (not pre_nan_.count(x)) eq % data[x];
            if (piracy_data and not piracy_nan_.count(x)) eq % (data["piracy"] * data[x]);
            if (public_data and not public_nan_.count(x)) eq % (data["public"] * data[x]);
        }
        // These don't get included in "pre":
        for (auto &x : {"param.cost_piracy", "param.piracy_link_proportion"}) {
            if (not data.has(x)) continue;
            if (piracy_data and not piracy_nan_.count(x)) eq % (data["piracy"] * data[x]);
            if (public_data and not public_nan_.count(x)) eq % (data["public"] * data[x]);
        }
        // These don't get included in "pre" or "piracy":
        for (auto &x : {"param.public_sharing_tax"}) {
            if (not data.has(x)) continue;
            if (public_data and not public_nan_.count(x)) eq % (data["public"] * data[x]);
        }
        marg_effects.add(eq);
    }

    marg_effects.solve();
    std::cout << "\n\n\nMarginal effects:\n=================\n" << marg_effects;
}

