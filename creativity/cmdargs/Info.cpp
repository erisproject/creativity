#include "creativity/cmdargs/Info.hpp"
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/errors.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <iostream>
#include <algorithm>
#include <string>
#include <cctype>

namespace creativity { namespace cmdargs {

namespace po = boost::program_options;

void Info::addOptions() {

    CmdArgs::addOptions(); // for --help, --version

    options_.add_options()
        ("precision,p", range<2,std::numeric_limits<double>::max_digits10>(output_precision), "      Specifies the precision level for result values.  The default is 6.")
        ("thin-periods,t", min<1>(thin_periods), "How many periods to show when --show-periods is specified.  1 means show every period, 2 means show every 2nd period, etc.  The default is 10.")
        ("all-periods,a", "Equivalent to --thin-periods=1, this shows all periods.")
        ("hide-periods,q", "If specified, hides the display of period summaries entirely.")
        ("show-cli-arguments,c", value(show_cli_args), "If specified, output the command line arguments to recreate the experiment.")
        ("memory-xz,M", value(memory_xz), "If the input file is an xz-compressed file, using this flag causes it to be decompressed into memory instead of writing it to a temporary file.")
        ("memory-all", value(memory), "If specified, copy the input file into memory before processing.  Implies --memory-xz")
        ;

    po::options_description input_desc("Input file");
    input_desc.add_options()
        ("input-file", value(input), "An input file as produced by creativity-cli or creativity-gui for which to display simulation details.")
        ;
    invisible_.add(input_desc);
    positional_.add("input-file", -1);
}

void Info::postParse(boost::program_options::variables_map &vars) {
    if (input.empty()) {
        throw po::required_option("FILENAME");
    }
    if (vars.count("all-periods")) thin_periods = 1;
    else if (vars.count("hide-periods")) thin_periods = 0;
}

std::string Info::usage() const {
    return CmdArgs::usage() + " FILENAME";
}

std::string Info::help() const {
    return CmdArgs::help() + "Input file:\n" +
        "  FILENAME                     A mandatory .crstate data file as produced by\n" +
        "                               creativity-cli or creativity-gui for which to\n" +
        "                               display simulation details.\n\n";
}

std::string Info::versionSuffix() const {
    return " -- simulation detail viewer";
}

}}
