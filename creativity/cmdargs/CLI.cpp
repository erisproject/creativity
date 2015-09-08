#include "creativity/cmdargs/CLI.hpp"
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/errors.hpp>
#include <boost/program_options/value_semantic.hpp>

namespace creativity { struct CreativitySettings; }

namespace creativity { namespace cmdargs {

namespace po = boost::program_options;

CLI::CLI(CreativitySettings &s) : Simulator(s) {}

// Single-letter options used:
// b B c C d D e f g G i j k K m M n N o O p P Q r R s S T u U w W x y z Z
// Available:
// a A E F H I J l L q t v V X Y

void CLI::addOptions() {
    output = "creativity-SEED.crstate";
    threads = 0;

    Simulator::addOptions();

    sim_controls_.add_options()
        ("output,o", value(output), "Output file for simulation results.  If this contains the characters 'SEED', they will be replaced with the random seed value used for the simulation.")
        ("tmpdir", value(tmpdir), "Output directory in which to write the output file while running the simulation.  When "
            "the simulation finishes, the temporary file is moved to the output location specified by -o.  If this argument is omitted, the "
            "file is written directly to its final location.")
        ("overwrite,O", "Allows output file given to -o to be overwritten.")
        ;
    options_.add(sim_controls_);
}

void CLI::postParse(boost::program_options::variables_map &vars) {
    Simulator::postParse(vars);
    overwrite = vars.count("overwrite") > 0;
}

}}
