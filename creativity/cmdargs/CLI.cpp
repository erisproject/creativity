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
        ("xz", value(xz), "Enable xz compression of simulation results file.  The filename specified to --output will have .xz appended.")
        ("memory-xz", value(memory_xz), "Implies --xz, but also keeps all results in memory until writing the final .xz file instead of writing to a temporary, uncompressed file first.  Warning: this option requires significantly more memory.")
        ("tmpdir", value(tmpdir), "Output directory in which to write the output file while running the simulation.  When "
            "the simulation finishes, the temporary file is moved to the output location specified by -o.  If this argument is omitted, the "
            "file is written to a temporary file in the same directory as the final output file.  Has no effect when --memory-xz is enabled.")
        ("overwrite,O", value(overwrite), "Allows output file given to -o to be overwritten.")
        ;
    options_.add(sim_controls_);
}

void CLI::postParse(boost::program_options::variables_map &vars) {
    Simulator::postParse(vars);
    if (memory_xz) xz = true;
}

std::string CLI::versionSuffix() const {
    return " -- command-line simulator";
}

}}
