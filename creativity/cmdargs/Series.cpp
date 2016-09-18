#include "creativity/cmdargs/Series.hpp"
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/errors.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <iostream>

namespace creativity { namespace cmdargs {

namespace po = boost::program_options;

Series::Series() {}

void Series::addOptions() {

    CmdArgs::addOptions(); // for --help, --version

    options_.add_options()
        ("series,s", value(series), "A comma-separated list of data series to calculate.  Run --help-series for a detailed list of possible options.")
        ("periods,T", value(periods), "    Specifies the number of periods for which variable series should be calculated.  The default, 0, uses all periods found in the first simulation file loaded.  All files must contain this number of periods.")
        ("allow-unused-periods,U", value(allow_unused_periods), "If specified, allow files to contain more than the number of periods specified with --periods.  If not given, an error occurs if files do not have matching period counts.  Has no effect if --periods is not explicitly specified.")
        ("piracy-begins,P", value(piracy_begins), "    Specifies the required period for piracy beginning.  The default, 0, obtains this from the first simulation file loaded.  All simulations must have the same value.")
        ("policy-begins,G", value(policy_begins), "Specifies the required period for policy beginning.  The default, 0, obtains this from the first simulation file loaded.  All simulations must have the same value.")
        ("ignore-errors,I", value(ignore_errors), "If specified, just warn about and skip any files that can't be read or are not compatible with the given or implied total, piracy, or policy periods.")
        ("precision", range<3,std::numeric_limits<double>::max_digits10>(double_precision), "      Specifies the precision level for floating point values.  The default is the minimum required to exactly represent all possible double values without any loss of precision.")
        ("output-directory,o", value(output_dir), "The directory in which to place series files.  If this directory does not exist, it will be created.  Each series S is written to a file named `series-S.csv' in this directory.  Existing files will be overwritten.")
        ("threads,j", value(threads), "    Maximum number of threads to use for data parsing.  0 (the default) disables data parsing threading entirely.")
        ("memory-xz,M", value(memory_xz), "If an input file is an xz-compressed file, using this flag causes it to be decompressed into memory instead of writing it to a temporary file.")
        ("tmpdir", value(tmpdir), "If --memory-xz is not specified, this specifies a temporary directory in which to place temporary decompressed files.  If omitted, the file is in the same directory as the input file.")
        ("help-series,S", value(help_series), "Shows the series that can be generated.")
        ;

    po::options_description input_desc("Input files");
    input_desc.add_options()
        ("input-files", value(input), "One or more input files from which to generate series.")
        ;
    invisible_.add(input_desc);
    positional_.add("input-files", -1);
}

void Series::postParse(boost::program_options::variables_map&) {
    if (output_dir.empty())
        throw std::logic_error("Error: --output-directory/-o argument is required");
}

std::string Series::usage() const {
    return CmdArgs::usage() + " FILE [FILE ...]";
}

std::string Series::help() const {
    return CmdArgs::help() + "Input files:\n" +
        "  FILE [FILE ...]                       One or more input files from which to\n" +
        "                                        calculate data series.  At least one must\n" +
        "                                        be given.\n\n" +
        "This program calculates the average values of various simulation variables in\n" +
        "each period, orders them across periods, then writes the results to CSV files.\n" +
        "The files can then be used to plot various quantiles of the simulation\n" +
        "parameters.\n\n";
}

std::string Series::versionSuffix() const {
    return " -- simulation data series generator";
}

}}
