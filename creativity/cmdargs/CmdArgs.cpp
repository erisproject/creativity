#include "creativity/cmdargs/CmdArgs.hpp"
#include "creativity/config.hpp"
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/errors.hpp>
#include <cstdlib>
#include <string>
#include <iostream>

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



void CmdArgs::parse(int argc, char const* const* argv) {
    if (options_.options().empty()) addOptions();

    prog_name_ = argv[0];

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
    version << "Creativity simulator v" << VERSION[0] << "." << VERSION[1] << "." << VERSION[2];
    return version.str();
}

std::string CmdArgs::usage() const {
    std::string name = prog_name_;
    if (name.empty()) name = "creativity";
    return "Usage: " + name + " [ARGS]";
}

std::string CmdArgs::help() const {
    std::ostringstream help;
    help << usage() << "\n";
    if (not options_.options().empty())
        help << "Supported arguments:\n" << options_ << "\n";
    return help.str();
}
}}
