#include "creativity/cmdargs/CLI.hpp"
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/errors.hpp>
#include <boost/program_options/value_semantic.hpp>

namespace creativity { struct CreativitySettings; }

namespace creativity { namespace cmdargs {

namespace po = boost::program_options;

CLI::CLI(CreativitySettings &s) : Simulator(s) {}

void CLI::addOptions() {
    output = "creativity-SEED.crstate";
    threads = 0;

    Simulator::addOptions();

    sim_controls_.add_options()
        ("quiet,q", value(quiet), "If specified, don't output current simulation status.")
        ("output,o", value(output), "Output file for simulation results.  If this contains the characters 'SEED', they will be replaced with the random seed value used for the simulation.")
        ("memory", value(memory), "Keeps all results in memory until writing the final crstate file instead of writing to an intermediate, uncompressed file first.  Warning: this option requires significantly more memory.")
        ("overwrite,O", value(overwrite), "Allows output file given to -o to be overwritten.")
        ;
    options_.add(sim_controls_);
}

void CLI::postParse(boost::program_options::variables_map &vars) {
    Simulator::postParse(vars);
}

std::string CLI::versionSuffix() const {
    return " -- command-line simulator";
}

}}
