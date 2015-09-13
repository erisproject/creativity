#include "creativity/cmdargs/Results.hpp"
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

void Results::addOptions() {

    CmdArgs::addOptions(); // for --help, --version

    options_.add_options()
        ("text,t", value(format_text), "Output the results in plain text.  This is the default when neither --html or --latex are specified.")
        ("html,h", value(format_html), "Output the results in HTML tables.  This cannot be used with --text or --latex.")
        ("latex,l", value(format_latex), "Output the results in LaTeX tables.  This cannot be used with --text or --html.")
        ("precision", range<2,std::numeric_limits<double>::max_digits10>(output_precision), "      Specifies the precision level for result values.  The default is 6.")
        ;

    po::options_description input_desc("Input file");
    input_desc.add_options()
        ("input-file", value(input), "An input file as produced by creativity-data from which to calculate estimators and related statistics.")
        ;
    invisible_.add(input_desc);
    positional_.add("input-file", -1);
}

void Results::postParse(boost::program_options::variables_map&) {
    if ((format_text and format_html) or (format_text and format_latex) or (format_html and format_latex)) {
        throw std::logic_error("Only one of --text, --html, and --latex may be specified.");
    }
    if (input.empty()) {
        throw po::required_option("FILENAME");
    }
    format = format_html ? Format::HTML :
        format_latex ? Format::LaTeX :
        Format::Text;
}

std::string Results::usage() const {
    return CmdArgs::usage() + " FILENAME";
}

std::string Results::help() const {
    return CmdArgs::help() + "Input file:\n" +
        "  FILENAME                     A mandatory data file as produced by\n" +
        "                               creativity-data for which to analyze results.\n\n";
}

std::string Results::versionSuffix() const {
    return " -- results analyzer";
}

}}
