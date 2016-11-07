#include "creativity/data/CSVParser.hpp"
#include "creativity/cmdargs/SeriesGraphs.hpp"
#include "creativity/data/quantiles.hpp"
#include "creativity/data/graph/Series.hpp"
#include "creativity/data/graph/PDF.hpp"
#include <eris/types.hpp>
#include <boost/filesystem/operations.hpp>
#include <cerrno>
#include <exception>
#include <string>
#include <iostream>
#include <regex>
#include <algorithm>

namespace creativity { namespace state { class State; } }

using namespace creativity;
using namespace creativity::data;
using namespace creativity::data::graph;
using namespace eris;

namespace fs = boost::filesystem;

// Struct for storing data parsing results
struct quantile_data {
    std::string filename;
    std::string title;
    // Map level -> { t -> [lower, upper] }, with levels in descending order
    std::map<double, std::map<int, std::pair<double, double>>, std::greater<double>> quantiles;
    // The median (i.e. the 0.5 quantile); map contains {t -> median} pairs.
    std::map<int, double> median;
    // Number of observations with finite data we extracted
    unsigned nobs = 0;
    // Used to track largest values found in this file:
    unsigned tmin = std::numeric_limits<unsigned>::max(),
             tmax = std::numeric_limits<unsigned>::min();
    double ymin = std::numeric_limits<double>::quiet_NaN(),
           ymax = std::numeric_limits<double>::quiet_NaN();

    void update_y_min_max(double min, double max) {
        if (std::isnan(ymin) or min < ymin) ymin = min;
        if (std::isnan(ymax) or max > ymax) ymax = max;
    }

    quantile_data(std::string file, std::string title) : filename{std::move(file)}, title{std::move(title)} {}
};

void extract_quantiles(CSVParser &parser, std::set<double> levels, quantile_data &in, int min_t, int max_t, bool quantiles_file) {
    std::unordered_map<double, size_t> l_pos_low, l_pos_high; // Quantile positions (in quantiles file)

    if (quantiles_file) {
        // Make sure we have all the quantiles we need, and store their row indices
        for (const auto &l : levels) {
            double qlow = 0.5 - l/2, qhigh = 0.5 + l/2;
            for (double q : {qlow, qhigh}) {
                auto qheader = quantile_field(q);
                if (not parser.hasField(qheader))
                    throw std::invalid_argument("file is missing required " + qheader + " quantile field");
            }
            l_pos_low.emplace (l, parser.fieldPosition(quantile_field(qlow)));
            l_pos_high.emplace(l, parser.fieldPosition(quantile_field(qhigh)));
        }
    }

    for (auto &row : parser) {
        std::vector<double> sorted_series_values;
        if (not quantiles_file) {
            // Extract finite values (ignore any NaNs or infinities)
            std::copy_if(row.begin()+1, row.end(), std::back_inserter(sorted_series_values), [](const double &v) { return std::isfinite(v); });
            if (sorted_series_values.empty()) continue; // Completely skip time periods with zero finite values (often initial rows, possibly others)
            std::sort(sorted_series_values.begin(), sorted_series_values.end());
        }

        long t = row[0];
        if ((min_t == -1 or t >= min_t) and (max_t == -1 or t <= max_t)) {
            if (t < in.tmin) in.tmin = t;
            if (t > in.tmax) in.tmax = t;
            in.nobs++;
            for (const auto &lvl : levels) {
                double low, high;
                if (quantiles_file) {
                    low = row[l_pos_low[lvl]];
                    high = row[l_pos_high[lvl]];
                }
                else {
                    low = creativity::data::quantile(sorted_series_values, 0.5 - lvl/2);
                    high = creativity::data::quantile(sorted_series_values, 0.5 + lvl/2);
                }

                in.update_y_min_max(low, high);
                if (lvl == 0) // Median
                    in.median.emplace(int(t), std::move(low));
                else // CI boundary pair
                    in.quantiles[lvl].emplace(int(t), std::make_pair(std::move(low), std::move(high)));
            }
        }
    }
}

void extract_series_ranges(CSVParser &parser, std::set<double> levels, quantile_data &in, int min_t, int max_t) {
    unsigned nfiles = parser.fields().size()-1;
    std::vector<double> score(nfiles, 0.);
    std::vector<std::vector<double>> rows;
    size_t nrow = 0;
    for (auto &row : parser) {
        long t = nrow++;
        if (row[0] != t) throw std::invalid_argument("series file has non-sequential t=" + std::to_string(row[0]) + " value on line " + std::to_string(parser.lineNumber()));

        if ((min_t != -1 and t < min_t) or (max_t != -1 and t > max_t)) {
            rows.emplace_back();
            continue; // Don't want this row
        }

        // Store the original row values so we can use them later
        rows.emplace_back(row.begin()+1, row.end());
        const auto &values = rows.back();
        // Sanity check:
        if (values.size() != nfiles)
            throw std::invalid_argument("series file has incorrect number of values on line " + std::to_string(parser.lineNumber()));

        // Copy the finite elements of the row, and sort them
        std::vector<double> sorted;
        std::copy_if(values.begin(), values.end(), std::back_inserter(sorted), [](const double &v) { return std::isfinite(v); });
        if (sorted.empty()) continue; // Everything gets the same (worst) score for this t, so don't bother adding it.
        std::sort(sorted.begin(), sorted.end());

        // Only count this as a useful observation if it has at least one finite value
        if (t < in.tmin) in.tmin = t;
        if (t > in.tmax) in.tmax = t;
        in.nobs++;

        // Calculate and add the marginal score for each source file
        double base = 0.5*(nfiles - sorted.size());
        for (size_t i = 0; i < values.size(); i++) {
            auto &v = values[i];
            auto &s = score[i];
            if (std::isfinite(v)) {
                // Find the first location in the sorted list that is bigger than v:
                auto found_it = std::find_if(sorted.begin(), sorted.end(), [&v](const double &check) { return check > v; });
                s += (base + std::distance(sorted.begin(), found_it) - 1) / (double) nfiles;
            }
            else {
                s += 0.25; // Maximum badness
            }
        }
    }

    auto sorted_scores = score;
    std::sort(sorted_scores.begin(), sorted_scores.end());
    for (const auto &l : levels) {
        double max_score = sorted_scores[std::lround(l*nfiles)];
        for (size_t t = in.tmin; t <= in.tmax; t++) {
            double max = std::numeric_limits<double>::quiet_NaN(), min = max;
            for (size_t i = 0; i < nfiles; i++) {
                double &v = rows[t][i];
                if (score[i] <= max_score and std::isfinite(v)) {
                    if (std::isnan(max) or v > max) max = v;
                    if (std::isnan(min) or v < min) min = v;
                    if (l == 0) {
                        // want median, and max_score is the smallest score, so we're done.  There
                        // might be duplicates, but that's okay: we're going with the first one we
                        // find.
                        break;
                    }
                }
            }

            if (std::isfinite(min)) {
                in.update_y_min_max(min, max);
                if (l == 0)
                    in.median.emplace(int(t), std::move(min));
                else
                    in.quantiles[l].emplace(int(t), std::make_pair(std::move(min), std::move(max)));
            }
        }
    }
}

int main(int argc, char *argv[]) {
    cmdargs::SeriesGraphs args;
    try {
        args.parse(argc, argv);
    }
    catch (const std::exception &e) {
        std::cerr << "\n" << e.what() << "\n\n";
        exit(5);
    }

    std::istringstream iss(args.levels);
    std::set<double> levels;
    bool invalid = false;
    std::string level;
    while (std::getline(iss, level, ',')) {
        std::size_t pos;
        if (level == "median") levels.insert(0);
        else {
            try {
                double q = std::stod(level, &pos);
                if (pos != level.size() or q < 0 or q > 1) throw std::invalid_argument("");
                levels.insert(q);
            } catch (const std::invalid_argument&) {
                std::cerr << "Error: requested level `" << level << "' is invalid\n";
                invalid = true;
            }
        }
    }

    if (invalid) {
        std::cerr << "Invalid confidence level(s) provided; aborting.\n\n";
        exit(1);
    }
    else if (levels.empty()) {
        std::cerr << "No confidence levels provided; aborting.\n\n";
        exit(1);
    }

    // Here's what we do now:
    // - Load each input file
    // - For series, calculate the needed quantiles
    // - Using the largest confidence level, calculate the minimum (from the lower quantile) and maximum (from the upper quantile).
    // - Find the lowest and highest "t" observations with data
    // - Repeat the above two steps across all input files
    // - This gives the x region (t range) and y region (data range) for the plots.
    // - For each file, plot the data:
    //   - In descending order of confidence levels, plot the ribbon associated with the confidence
    //     level; for 0 (the median) just plot a line instead.

    std::list<quantile_data> input;

    // DEBUG:
    auto trim = common_ends(args.input.begin(), args.input.end());

    for (const auto &file : args.input) {
        fs::path input_path(file);
        if (not fs::is_regular_file(input_path)) {
            std::cerr << "Error: input file `" << file << "' does not exist or is not a regular file\n";
            exit(3);
        }

        std::ostringstream title;
        for (auto it = args.title.begin(); it != args.title.end(); it++) {
            if (*it == '%') {
                if (++it != args.title.end()) {
                    switch (*it) {
                        case '%': title << '%'; break;
                        case 'F': title << file; break;
                        case 'p': title << input_path.parent_path(); break;
                        case 'f': title << input_path.filename(); break;
                        case 'u': title << file.substr(trim.first, file.length()-trim.second-trim.first); break;
                        case 'w': title << std::regex_replace(input_path.filename().replace_extension().string(), std::regex(".*\\W"), ""); break;
                        default: title << '%' << *it; // Not actually an escape; pass both % and the value through
                    }
                }
            }
            else title << *it;
        }

        CSVParser parser(file);
        input.emplace_back(file, title.str());
        auto &in = input.back();

        int min_t = std::max(args.data_min_t, args.graph_min_t);
        int max_t = args.data_max_t == -1 ? args.graph_max_t : args.graph_max_t == -1 ? args.data_max_t :
            std::min(args.data_max_t, args.graph_max_t);

        try {
            if (parser.fields()[0] != "t") throw std::invalid_argument("first field != t");

            // First figure out whether this is a series or quantiles file:
            bool found_series = false, found_quantiles = false;
            for (size_t i = 1; i < parser.fields().size(); i++) {
                if (std::regex_match(parser.fields()[i], quantile_field_regex))
                    found_quantiles = true;
                else
                    found_series = true;

                // Both isn't allowed!
                if (found_quantiles and found_series)
                    throw std::invalid_argument("file headers invalid: found both quantile and non-quantile headers");
            }
            if (not found_series and not found_quantiles)
                throw std::invalid_argument("file headers do not appear to be either a series or quantiles file");

            if (found_quantiles and not args.per_t_confidence)
                throw std::invalid_argument("Can't use --source-confidence with quantiles file");

            if (args.per_t_confidence) {
                extract_quantiles(parser, levels, in, min_t, max_t, found_quantiles);
            }
            else {
                extract_series_ranges(parser, levels, in, min_t, max_t);
            }

            // Warn if the admitted observations seems too low
            if (in.nobs < 10) {
                std::cerr << "Warning: input file `" << file << "' had " << in.nobs << " admissable observations (perhaps --t-{min,max,from,to} values are too restrictive?)";
            }
        }
        catch (const std::invalid_argument &e) {
            std::cerr << "\n\nError: input file `" << file << "' is not a valid series or quantiles file: " << e.what() << ".  Aborting.\n\n";
            exit(8);
        }
    }

    unsigned tmin = std::numeric_limits<unsigned>::max(), tmax = std::numeric_limits<unsigned>::min();
    double ymin = std::numeric_limits<double>::infinity(), ymax = -std::numeric_limits<double>::infinity();
    for (const auto &f : input) {
        if (f.tmin < tmin) tmin = f.tmin;
        if (f.tmax > tmax) tmax = f.tmax;
        if (f.ymin < ymin) ymin = f.ymin;
        if (f.ymax > ymax) ymax = f.ymax;
    }

    if (args.graph_min_t != -1) tmin = args.graph_min_t;
    if (args.graph_max_t != -1) tmax = args.graph_max_t;
    if (tmin >= tmax) throw std::logic_error("Error: minimum t >= maximum t; cannot continue with graph.  Check --t-{min,max,to,from} arguments.");
    if (std::isfinite(args.graph_min_value)) ymin = args.graph_min_value;
    if (std::isfinite(args.graph_max_value)) ymax = args.graph_max_value;
    if (ymin >= ymax) throw std::logic_error("Error: minimum data value >= maximum data value; cannot continue with graph.  Check --y-{min,max} arguments.");

    fs::path output_path(args.output);
    if (output_path.extension() != ".pdf") {
        // FIXME
        std::cerr << "Support for non-PDF output files not yet implemented";
        exit(10);
    }
    if (not args.overwrite and fs::exists(output_path)) {
        std::cerr << "Error: output file " << output_path << " already exists; specify another output file, or use --overwrite to overwrite\n";
        exit(11);
    }

    PDF output(args.output, args.width, args.height);
    Series graph(output, tmin, tmax, ymin, ymax);
    graph.tick_grid_style.colour = RGBA(0.9);
    graph.tick_grid_extra_style.colour = RGBA(0.5);
    graph.legend_position = args.legend_position;
    graph.legend_rel_x = args.legend_rel_x;
    graph.legend_rel_y = args.legend_rel_y;
    graph.x_label = "t";
    FillStyle confband_style = Series::default_region_style;
    LineStyle median_style = Series::default_line_style;
    {
        // Figure out the legend by multiplying colours together
        auto legend_style = Series::default_legend_box_style;
        legend_style.fill_colour = graph.graph_style.fill_colour;
        RGBA bottom = legend_style.fill_colour * confband_style.fill_colour;
        for (auto it = levels.rbegin(); it != levels.rend(); it++) {
            if (*it > 0) {
                std::ostringstream ci;
                ci << (*it * 100) << "% conf. boundary";
                graph.addLegendItem(ci.str(), bottom, true, legend_style);
                legend_style.fill_colour = bottom;
                bottom *= confband_style.fill_colour;
            }
            else {
                graph.addLegendItem("Median", median_style, true, legend_style);
            }
        }
    }

    bool first = true;
    for (const auto &in : input) {
        if (first) first = false;
        else graph.newPage();
        int local_tmin = tmin, local_tmax = tmax;
        double local_ymin = ymin, local_ymax = ymax;
        if (not args.same_horizontal_scale) {
            local_tmin = args.graph_min_t >= 0 ? args.graph_min_t : in.tmin;
            local_tmax = args.graph_max_t >= 0 ? args.graph_max_t : in.tmax;
        }
        if (not args.same_vertical_scale) {
            local_ymin = std::isfinite(args.graph_min_value) ? args.graph_min_value : in.ymin;
            local_ymax = std::isfinite(args.graph_max_value) ? args.graph_max_value : in.ymax;
        }
        graph.setExtents(local_tmin, local_tmax, local_ymin, local_ymax, false);
        graph.recalcTicks(10, 8, Series::TickEnds::Add);
        graph.title = in.title;
        graph.t_grid_extra.insert(args.t_extra.cbegin(), args.t_extra.cend());

        for (const auto &conf : in.quantiles)
            graph.addRegion(conf.second, confband_style);
        if (not in.median.empty())
            graph.addLine(in.median, median_style);
    }
}

