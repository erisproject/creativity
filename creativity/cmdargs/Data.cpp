#include "creativity/cmdargs/Data.hpp"
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/errors.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <iostream>

namespace creativity { namespace cmdargs {

namespace po = boost::program_options;

Data::Data() {}

void Data::addOptions() {

    CmdArgs::addOptions(); // for --help, --version

    options_.add_options()
        ("human-readable,H", value(human_readable), "Produce output in human-readable format.  The default (when this is not specified) outputs comma-separated values.  This also changes the default value for --precision to 8.")
        ("no-csv-header", value(no_csv_header), "Don't include the CSV header.  Useful when combining output from multiple creativity-data invocations.")
        ("only-csv-header", value(only_csv_header), "Outputs just the CSV header associated with the program arguments instead of processing files.  This is primarily intended to be combined with other invocations using the --no-csv-header argument.")
        ("precision", range<3,std::numeric_limits<double>::max_digits10>(double_precision), "      Specifies the precision level for floating point values.  The default is the minimum required to exactly represent all possible double values without any loss of precision.")
        ("periods,t", min<1>(data_periods), "    Specifies the number of periods to use for calculating pre, new, and post-piracy data")
        ("skip-piracy", value(skip.piracy), "If specified, do not produce data for piracy periods.  Required for simulation data that does not contain pre-public sharing piracy periods.")
        ("skip-public-sharing", value(skip.public_sharing), "If specified, do not produce data for public sharing periods.  Required for simulation data that does not contain public sharing periods.")
        ("skip-short-run", value(skip.short_run), "If specified, don't include \"short-run\" effects generated from the initial periods of piracy or public sharing.")
        ("verify-periods,T", value(verify.periods), "    If non-zero, only use given simulation files with the last time period matching the given argument.")
        ("verify-piracy-begins,P", value(verify.piracy_begins), "If non-zero, only use given simulation files with piracy beginning in the given period.")
        ("verify-public-sharing-begins,G", value(verify.public_sharing_begins), "If non-zero, only use given simulation files with public sharing beginning in the given period.")
        ("threads,j", value(threads), "    Maximum number of threads to use for data parsing.  0 (the default) disables data parsing threading entirely.")
        ("memory-xz,M", value(memory_xz), "If an input file is an xz-compressed file, using this flag causes it to be decompressed into memory instead of writing it to a temporary file.")
        ("tmpdir", value(tmpdir), "If --memory-xz is not specified, this specifies a temporary directory in which to place the temporary decompressed file.  If omitted, the file is in the same directory as the input file.")
        ;

    po::options_description input_desc("Input files");
    input_desc.add_options()
        ("input-files", value(input), "One or more input files from which to calculate data.")
        ;
    invisible_.add(input_desc);
    positional_.add("input-files", -1);
}

void Data::postParse(boost::program_options::variables_map &vars) {
    if (human_readable and vars["precision"].defaulted()) {
        double_precision = 8;
    }
    if (only_csv_header and no_csv_header)
        throw std::logic_error("Conflicting arguments: --only-csv-header cannot be used with --no-csv-header");
    if (human_readable and (only_csv_header or no_csv_header))
        throw std::logic_error("Conflicting arguments: --human-readable cannot be used with --{no,only}-csv-header");
    if (only_csv_header and not input.empty())
        throw std::logic_error("Conflicting arguments: --only-csv-header cannot be used when input files are specified");
}

std::string Data::usage() const {
    return CmdArgs::usage() + " FILE [FILE ...]";
}

std::string Data::help() const {
    return CmdArgs::help() + "Input files:\n" +
        "  FILE [FILE ...]                       One or more input files from which to\n" +
        "                                        calculate statistics.  At least one must\n" +
        "                                        be given.\n\n";
}

std::string Data::versionSuffix() const {
    return " -- simulation data collector";
}

}}
