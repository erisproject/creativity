#include "creativity/cmdargs/SeriesGraphs.hpp"
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/errors.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <iostream>

namespace creativity { namespace cmdargs {

namespace po = boost::program_options;

SeriesGraphs::SeriesGraphs() {}

void SeriesGraphs::addOptions() {

    CmdArgs::addOptions(); // for --help, --version

    options_.add_options()
        ("levels,l", value(levels), "A comma-separated list of confidence levels to plot.  0 corresponds to just plotting the median, 1 corresponds to a range that includes all observations.  The default plots the median, inter-quartile region (.5), and the 90 and 95\% confidence intervals (.9, .95).")
        ("output,o", value(output), "The filename to output results.  Should be either FILENAME.pdf or FILENAME.png; for the latter, the first series will be output to FILENAME.png; subsequent series will be plotted in FILENAME-2.png, FILENAME-3.png, etc.")
        ("overwrite,O", value(overwrite), "If an output filename (after applying --prefix and --unprefix options) already exists, overwrite it.  Without this option, an error results if an output file already exists.")
        ("width,W", above<0>(width), "The width of the output file(s), in inches.")
        ("height,H", above<0>(height), "The height of the output file(s), in inches.")
        ("resolution,R", above<0>(resolution), "The resolution of the images for PNG output, in pixels per inch.")
        ;

    po::options_description input_desc("Input files");
    input_desc.add_options()
        ("input-files", value(input), "One or more series input files from which to generate a quantile output graphs.")
        ;
    invisible_.add(input_desc);
    positional_.add("input-files", -1);
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
