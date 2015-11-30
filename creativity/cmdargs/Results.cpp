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

    po::options_description output_desc("Output"), analysis_desc("Results analysis"), input_desc("Input file"), format_desc("Format options");

    analysis_desc.add_options()
        ("summary,S", value(analysis.summary), "Show data summary: number of simulations in total, and number in each category.")
        ("write-vs-nowrite,W", value(analysis.write_or_not), "Show parameter distributions conditional on simulations with consistent writing and simulations with no-writing stages.")
        ("write-vs-nowrite-corr,C", value(analysis.write_or_not_corrcov), "Include parameter correlations (only has effect when --write-vs-nowrite is activated); disabled by default.")
        ("average-effects,A", value(analysis.average), "Show average model effects across different simulation stages.")
        ("marginal-effects,M", value(analysis.marginal), "Show marginal model effects across different simulation stages and model parameters.")
        ("all,a", value(analysis.all), "Implies --summary --write-vs-nowrite --average-effects --marginal-effects, but not --write-vs-nowrite-corr.  If none of the above nor --none are given, this is the default.")
//        ("none,n", value(analysis.none), "Show none of the above analysis.  May not be combined with any of this above.  This flag is intended for use with the --dump-* options.")
        ;
    options_.add(analysis_desc);

    format_desc.add_options()
        ("text,t", value(format.text), "Output the results in plain text.  This is the default when neither --html or --latex are specified.")
        ("html,h", value(format.html), "Output the results in HTML tables.  This cannot be used with --text or --latex.")
        ("latex,l", value(format.latex), "Output the results in LaTeX tables.  This cannot be used with --text or --html.")
        ("precision", range<2,std::numeric_limits<double>::max_digits10>(format.precision), "      Specifies the precision level for result values.  The default is 6.")
        ;
    options_.add(format_desc);

    output_desc.add_options()
        ("output,o", value(output.filename), "If specified, data analysis results will be written to the given filename instead of displayed.  The file must not exist unless --overwrite is specified.")
        //("dump-write", value(output.dump_writing), "If specified, dump intermediate data of simulations with writing in each stage to the given CSV file.  The file must not exist unless --overwrite is specified.")
        //("dump-nowrite", value(output.dump_nowriting), "If specified, dump intermediate data of simulations *without* writing in one or more stages to the given CSV file.  The file must not exist unless --overwrite is specified.")
        //("overwrite,O", value(output.overwrite), "If specified, files given to --output, --dump-write or --dump-nowrite will be overwritten if they exist.")
        ("overwrite,O", value(output.overwrite), "If specified, filename given to --output will be overwritten if it exists.")
        ("no-preamble,P", value(output.no_preamble), "Suppress preamble/postamble for HTML or LaTeX output formats.  The default, without this option, includes preamble/postamble output that makes the output a valid HTML/LaTeX document; with this option, only the relevent document fragment is output.");
        ;
    options_.add(output_desc);

    input_desc.add_options()
        ("input-file", value(input), "An input file as produced by creativity-data from which to calculate estimators and related statistics.")
        ;
    invisible_.add(input_desc);
    positional_.add("input-file", -1);
}

void Results::postParse(boost::program_options::variables_map&) {
    if (format.text + format.html + format.latex > 1) {
        throw std::logic_error("Only one of --text, --html, and --latex may be specified.");
    }
    format.type = format.html ? data::TableFormat::HTML :
        format.latex ? data::TableFormat::LaTeX :
        data::TableFormat::Text;

    if (input.empty()) {
        throw po::required_option("FILENAME");
    }

    if (analysis.none) {
        if (analysis.all or analysis.summary or analysis.write_or_not or analysis.average or analysis.marginal)
            throw std::logic_error("Option --none cannot be combined with --all, --write-vs-nowrite, --average-effects, --marginal-effects");
    }
    else if (
            analysis.all or // --all explicitly given
            (not analysis.write_or_not and not analysis.average and not analysis.marginal)) // Nothing given: --all is default
        analysis.summary = analysis.write_or_not = analysis.average = analysis.marginal = true;
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
