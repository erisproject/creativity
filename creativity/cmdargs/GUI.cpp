#include "creativity/cmdargs/GUI.hpp"
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/errors.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <thread>

namespace creativity { namespace cmdargs {

namespace po = boost::program_options;

GUI::GUI(CreativitySettings &s) : Simulator(s) {}

// Single-letter options used (bracketed values are in CLI but not GUI):
// b B c C d D e f g G i j k K m M n N o [O] p P Q r R s S T u U w W x y z Z
// Available:
// a A E F H I J l L q t v V X Y

void GUI::addOptions() {
    threads = std::thread::hardware_concurrency();
    if (threads == 1) threads = 0;

    Simulator::addOptions();

    sim_controls_.add_options()
        ("output,o", value(output), "Output file for simulation results.  If this contains the characters 'SEED', they will be replaced with the random seed value used for the simulation.  If omitted, the simulation is not written to disk.")
        ("initialize", "If specified, initialize the simulation using the given settings, but do no start it.  Otherwise the GUI starts uninitialized, but ready-to-run with the given settings.  Ignored if --start is specified.")
        ("start", "If specified, start running immediately using the given settings.  Otherwise the GUI starts uninitialized, but ready-to-run with the given settings.")
        ;
    options_.add(sim_controls_);

    // FIXME: could add colour settings here

    po::options_description input_desc("Input file");
    input_desc.add_options()
        ("input-file", value(input), "Previous simulation to load instead of configuring a new simulation.")
        ;
    invisible_.add(input_desc);
    positional_.add("input-file", -1);
}


void GUI::postParse(boost::program_options::variables_map &vars) {
    Simulator::postParse(vars);
    start = vars.count("start") > 0;
    initialize = vars.count("initialize") > 0;
}

std::string GUI::usage() const {
    return Simulator::usage() + " [FILE-TO-LOAD]";
}

std::string GUI::help() const {
    return Simulator::help() + "Loading an existing simulation:\n" +
        "  FILE-TO-LOAD                          An existing crstate file to display\n" +
        "                                        instead of configuring a new simulation.\n\n";
}

std::string GUI::versionSuffix() const {
    return " -- graphical simulator and simulation viewer";
}

}}
