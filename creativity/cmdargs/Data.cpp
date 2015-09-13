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
        ("precision", range<3,std::numeric_limits<double>::max_digits10>(double_precision), "      Specifies the precision level for floating point values.  The default is the minimum required to exactly represent all possible double values without any loss of precision.")
        ("periods,t", min<1>(data_periods), "    Specifies the number of periods to use for calculating pre, new, and post-piracy data")
        ("skip-piracy", value(skip.piracy), "If specified, do not produce data for piracy periods.  Required for simulation data that does not contain pre-public sharing piracy periods.")
        ("skip-public-sharing", value(skip.public_sharing), "If specified, do not produce data for public sharing periods.  Required for simulation data that does not contain public sharing periods.")
        ("verify-periods,T", value(verify.periods), "    If non-zero, only use given simulation files with the last time period matching the given argument.")
        ("verify-piracy-begins,P", value(verify.piracy_begins), "If non-zero, only use given simulation files with piracy beginning in the given period.")
        ("verify-public-sharing-begins,G", value(verify.public_sharing_begins), "If non-zero, only use given simulation files with public sharing beginning in the given period.")
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
