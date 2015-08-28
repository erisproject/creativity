#include "creativity/cmdargs/CmdArgs.hpp"
#include "creativity/config.hpp"
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <eris/types.hpp>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <sstream>
#include <iostream>
#include <thread>

namespace po = boost::program_options;

namespace creativity { namespace cmdargs {

void CmdArgs::addOptions() {
    po::options_description help("About");

    help.add_options()
        ("help,h", "Displays usage information.")
        ("version", "Displays version information.")
        ;

    options_.add(help);
}

void CmdArgs::postParse(boost::program_options::variables_map&) {}



void CmdArgs::parse(int argc, char *argv[]) {
    if (options_.options().empty()) addOptions();

    po::options_description all_opts;
    all_opts.add(options_);
    all_opts.add(invisible_);

    // The variable map storing specified options; we don't keep this: everything gets set directly
    boost::program_options::variables_map vars;
    store(po::command_line_parser(argc, argv).options(all_opts).positional(positional_).run(), vars);
    notify(vars);

    bool need_help = vars.count("help") > 0, need_version = vars.count("version") > 0;
    if (need_help) {
        std::cout << version() << "\n\n" << help();
        std::exit(0);
    }
    else if (need_version) {
        std::cout << version() << "\n\n";
        std::exit(0);
    }

    postParse(vars);

}

std::string CmdArgs::version() const {
    std::ostringstream version;
    version << "Creativity simulator v" << VERSION[0] << "." << VERSION[1] << "." + VERSION[2];
    return version.str();
}

std::string CmdArgs::help() const {
    std::ostringstream help;
    help << options_ << "\n";
    return help.str();
}
}}
