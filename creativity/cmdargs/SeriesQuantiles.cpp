#include "creativity/cmdargs/SeriesQuantiles.hpp"
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/errors.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <iostream>

namespace creativity { namespace cmdargs {

namespace po = boost::program_options;

SeriesQuantiles::SeriesQuantiles() {}

void SeriesQuantiles::addOptions() {

    CmdArgs::addOptions(); // for --help, --version

    options_.add_options()
        ("quantiles,q", value(quantiles), "A comma-separated list of quantiles to calculate.  The default is 0,.005,.01,.025,.05,.25,.5,.75,.95,.975,.99,.995,1")
        ("precision", range<3,std::numeric_limits<double>::max_digits10>(double_precision), "      Specifies the precision level for floating point values.  The default is the minimum required to exactly represent all possible double values without any loss of precision.")
        ("prefix,p", value(output_prefix), "A prefix to prepend to input filenames to construct an output filename.")
        ("unprefix,P", value(output_unprefix), "If the given prefix is found in input filenames, it is removed before applying the --prefix value.")
        ("overwrite,O", value(overwrite), "If an output filename (after applying --prefix and --unprefix options) already exists, overwrite it.  Without this option, an error results if an output file already exists.")
        ;

    po::options_description input_desc("Input files");
    input_desc.add_options()
        ("input-files", value(input), "One or more series input files from which to generate quantile output files.")
        ;
    invisible_.add(input_desc);
    positional_.add("input-files", -1);
}

std::string SeriesQuantiles::usage() const {
    return CmdArgs::usage() + " FILE [FILE ...]";
}

std::string SeriesQuantiles::help() const {
    return CmdArgs::help() + "Input files:\n" +
        "  FILE [FILE ...]                       One or more series input files from which to\n" +
        "                                        calculate quantiles.  At least one must\n" +
        "                                        be given.\n\n" +
        "This program takes files as produced by the creativity-series program as inputs\n" +
        "and retains the existing sample quantiles for each period in the series files.\n\n";
}

std::string SeriesQuantiles::versionSuffix() const {
    return " -- simulation data series quantile converter";
}

}}
