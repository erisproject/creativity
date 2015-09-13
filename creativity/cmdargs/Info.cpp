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
        ("precision", range<2,std::numeric_limits<double>::max_digits10>(output_precision), "      Specifies the precision level for result values.  The default is 6.")
        ;

    po::options_description input_desc("Input file");
    input_desc.add_options()
        ("input-file", value(input), "An input file as produced by creativity-cli or creativity-gui for which to display simulation details.")
        ;
    invisible_.add(input_desc);
    positional_.add("input-file", -1);
}

void Info::postParse(boost::program_options::variables_map&) {
    if (input.empty()) {
        throw po::required_option("FILENAME");
    }
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
