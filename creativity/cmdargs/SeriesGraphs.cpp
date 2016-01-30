#include "creativity/cmdargs/SeriesGraphs.hpp"
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/errors.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <unordered_map>
#include <iostream>

namespace creativity { namespace cmdargs {

namespace po = boost::program_options;

SeriesGraphs::SeriesGraphs() {}

void SeriesGraphs::addOptions() {

    CmdArgs::addOptions(); // for --help, --version

    options_.add_options()
        ("levels,l", value(levels), "A comma-separated list of confidence levels to plot.  0 corresponds to the median, 1 corresponds to a range that includes all observations.  The default plots the median, inter-quartile region (.5), and the 90 and 95\% confidence intervals (.9, .95).")
        ("output,o", value(output), "The filename to output results.  Should be either FILENAME.pdf or FILENAME.png; for the latter, the first series will be output to FILENAME.png; subsequent series will be plotted in FILENAME-2.png, FILENAME-3.png, etc.")
        ("overwrite,O", value(overwrite), "If an output filename (after applying --prefix and --unprefix options) already exists, overwrite it.  Without this option, an error results if an output file already exists.")
        ("title,t", value(title), "A title for the generated graph.  The title supports the following substitutes: %F - input path and filename; %p - input directory; %f - input filename; %u - unique portion of input path/filename (across all input files); %w - the last pre-extension [a-zA-Z0-9_] characters of the input file; %% - a literal % symbol.")
        ("width,W", above<0>(width), "The width of the output file(s), in inches.")
        ("height,H", above<0>(height), "The height of the output file(s), in inches.")
        ("resolution,r", above<0>(resolution), "The resolution of the images for PNG output, in pixels per inch.")
        ("t-min", min<-1>(graph_min_t), "The t value at which to start the graph.  If -1, the graph starts at the lowest t value found in the input file(s).")
        ("t-max", min<-1>(graph_max_t), "The t value at which to end the graph.  If -1, the graph ends at the largest t value found in the input file(s).")
        ("y-min", value(graph_min_value), "The data value at which the y-axis should begin.  If NaN (the default), the y-axis is sized to fit all of the required values.")
        ("y-max", value(graph_max_value), "The data value at which the y-axis should end.  If NaN (the default), the y-axis is sized to fit all of the required values.")
        ("t-from", min<-1>(data_min_t), "The t value at which to start including observations.  If -1, include all observations.  In contrast to --t-min, this controls when data is plotted, not the graph axis limits.")
        ("t-to", min<-1>(data_max_t), "The t value at which to stop including observations.  If -1, include all observations.  In contrast to --t-max, this controls when data is plotted, not the graph axis limits.")
        ("same-t-scale,T", "If specified, the time period scale for all graphs will be the same.  The default determines each graph's t scale separately.  This option is implied if both --t-min and --t-max are >= 0.")
        ("different-y-scale,Y", "If specified, the value scale (y-axis) will be determined to fit the data in each file separately.  The default uses the same scale for all graphs.  This option cannot be used when both --y-min and --y-max values are given.")
        ("legend-position,L", value(legend_position_input), "The general position of the graph legend, one of 'none', 'inside', 'right', 'left', 'top', 'bottom' where 'none' suppresses the legend, 'inside' puts it inside the graph, and the rest put it outside on the indicated side of the graph.")
        ("legend-x", range<0,1>(legend_rel_x), "Relative x position (for --legend-position={inside,top,bottom}): 0 is left-most, 1 is right-most.")
        ("legend-y", range<0,1>(legend_rel_y), "Relative y position (for --legend-position={inside,left,right}): 0 is top-most, 1 is bottom-most.")
        ;

    po::options_description input_desc("Input files");
    input_desc.add_options()
        ("input-files", value(input), "One or more series input files from which to generate a quantile output graphs.")
        ;
    invisible_.add(input_desc);
    positional_.add("input-files", -1);
}

using LP = creativity::data::graph::Series::LegendPosition;
const std::unordered_map<std::string, LP> posmap{
    {"none", LP::None},
    {"right", LP::Right},
    {"left", LP::Left},
    {"bottom", LP::Bottom},
    {"top", LP::Top},
    {"inside", LP::Inside}
};

void SeriesGraphs::postParse(boost::program_options::variables_map &vars) {
    if (graph_min_t >= 0 and graph_max_t >= 0 and graph_max_t <= graph_min_t)
        throw po::invalid_option_value("--t-max value must be larger than --t-min value");
    if (std::isfinite(graph_min_value) and std::isfinite(graph_max_value) and graph_max_value <= graph_min_value)
        throw po::invalid_option_value("--y-max value must be larger than --y-min value");
    if (vars.count("same-t-scale") or (graph_min_t >= 0 and graph_max_t >= 0))
        same_horizontal_scale = true;
    if (vars.count("different-y-scale")) {
        if (std::isfinite(graph_min_value) and std::isfinite(graph_max_value))
            throw po::invalid_option_value("Option --different-y-scale cannot be used when both --y-min and --y-max values are given");
        same_vertical_scale = false;
    }
    if (not legend_position_input.empty()) {
        try {
            legend_position = posmap.at(legend_position_input);
        } catch (const std::out_of_range&) {
            throw po::invalid_option_value("Invalid value `" + legend_position_input + "' specified for --legend-position");
        }
    }
}

std::string SeriesGraphs::usage() const {
    return CmdArgs::usage() + " FILE [FILE ...]";
}

std::string SeriesGraphs::help() const {
    return CmdArgs::help() + "Input files:\n" +
        "  FILE [FILE ...]                       One or more series input files from which to\n" +
        "                                        calculate quantiles.  At least one must\n" +
        "                                        be given.\n\n" +
        "This program takes files as produced by the creativity-series program as inputs\n" +
        "and retains the existing sample quantiles for each period in the series files.\n\n";
}

std::string SeriesGraphs::versionSuffix() const {
    return " -- simulation data series quantile converter";
}

}}
