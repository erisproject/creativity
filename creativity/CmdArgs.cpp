#include "creativity/CmdArgs.hpp"
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

namespace creativity {

template <>
std::string CmdArgs::output_string(double v) {
    return std::regex_replace(
            std::regex_replace(std::to_string(v),
                std::regex("(\\.\\d*?)0+$"),
                "$1"),
            std::regex("\\.$"),
            "");
}

void CmdArgs::addOptions() {
    po::options_description help("About");

    help.add_options()
        ("help,h", "Displays usage information.")
        ("version", "Displays version information.")
        ;

    options_.add(help);
}

void CmdArgs::postParse(boost::program_options::variables_map&) {}

void CmdArgs::Simulation::addOptions() {
    CmdArgs::addOptions();
    po::options_description structure("Structure"), initial("Initial Behaviour"),
        costs("Costs"), beliefs("Beliefs"), piracy("Piracy"), publicprov("Public Provisioning");

    structure.add_options()
        ("dimensions,D", min<1>(s_.dimensions), "Number of dimensions of the simulation")
        ("readers,r", min<1>(s_.readers), "Number of reader/author agents in the simulation")
        ("density,d", above<0>(density_), (u8"Reader density (in readers per unit^D, where D is the configured # of dimensions).  The default is the density required to have simulation boundaries at ±" + output_string(s_.boundary) + " in each dimension.").c_str())
        ("reader-step-sd,R", min<0>(s_.reader_step_sd), "Standard deviation of the inter-period random-direction reader movement distance ~ |N(0, R)|")
        ("book-distance-sd,B", min<0>(s_.book_distance_sd), "Standard deviation of book distance from author; distance ~ |N(0, B)|")
        ("book-quality-sd,Q", min<0>(s_.book_quality_sd), "Standard deviation of book perceived quality; perceived quality ~ N(q, Q), where q is the innate quality and Q is this value")
        ("reader-creation-shape,s", below<1>(s_.reader_creation_shape), u8"Shape parameter, β, of the creator effort function: q(ℓ) = (α/β)[(ℓ+1)^β - 1]")
        ("reader-creation-scale-min,z", min<0>(s_.reader_creation_scale_min), u8"Minimum support of scale parameter α ~ U[a,b] for the creator effort function: q(ℓ) = (α/β)[(ℓ+1)^β - 1]")
        ("reader-creation-scale-max,Z", min<0>(s_.reader_creation_scale_max), u8"Maximum support of scale parameter α ~ U[a,b] for the creator effort function: q(ℓ) = (α/β)[(ℓ+1)^β - 1]")
        ("creation-time,e", value(s_.creation_time), u8"The number of periods that elapse between the creation decision and the book being ready")
        ;
    options_.add(structure);

    initial.add_options()
        ("initial-prob-write,x", range<0, 1>(s_.initial.prob_write), "The probability of writing in initial periods")
        ("initial-quality-min,m", min<0>(s_.initial.q_min), "The minimum support of quality q ~ U[a,b] for authored books in initial periods")
        ("initial-quality-max,M", min<0>(s_.initial.q_max), "The maximum support of quality q ~ U[a,b] for authored books in initial periods")
        ("initial-price-min,n", min<0>(s_.initial.p_min), "The minimum support of price c + U[a,b] for new books in initial periods")
        ("initial-price-max,N", min<0>(s_.initial.p_max), "The maximum support of price c + U[a,b] for new books in initial periods")
        ("initial-prob-keep,k", range<0, 1>(s_.initial.prob_keep), "The probability of keeping a previously-written book on the market for another period")
        ("initial-keep-price,K", range<0, 1>(s_.initial.keep_price), "The  price-above-marginal-cost level (relative to current P-MC) for a book left on the market for another period")
        ;
    options_.add(initial);

    costs.add_options()
        ("income,i", above<0>(s_.income), "Per-period external reader income")
        ("cost-fixed,C", min<0>(s_.cost_fixed), "Fixed cost of keeping a book on the market for a period")
        ("cost-unit,c", min<0>(s_.cost_unit), "Unit cost of making a copy of a book")
        ("cost-piracy,y", min<0>(s_.cost_piracy), "Cost of receiving a pirated copy of a book")
        ;
    options_.add(costs);

    beliefs.add_options()
        ("prior-scale,w", min<1>(s_.prior_scale), "The per-period standard deviation scaling factor applied when a previous belief is used as the next period's prior")
        ("prediction-draws,p", min<1>(s_.prediction_draws), "The number of draws agents take from beliefs when using for prediction.  Lower values can be used to introduce deliberate randomness into the simulation")
        ("burnin-periods,u", value(s_.burnin_periods), "The number of initial periods during which `--prior-scale-burnin' should be used instead of `--prior-scale'")
        ("belief-threshold,b", value(s_.initial.belief_threshold), "The minimum n-k value at which a readers bases decision on beliefs instead of initial parameters")
        ("prior-scale-burnin,U", min<1>(s_.prior_scale_burnin), "The same as --prior-weight, but applied in the first `--burnin-periods' periods")
        ;
    options_.add(beliefs);

    piracy.add_options()
        ("piracy-begins,P", value(s_.piracy_begins), "The period in which piracy becomes available")
        ("piracy-link-proportion,f", range<0, 1>(s_.piracy_link_proportion), "Proportion of potential sharing links between readers that are created")
        ("prior-scale-piracy,W", min<1>(s_.prior_scale_piracy), "The same as --prior-weight, but applied in the first piracy period")
        ;
    options_.add(piracy);

    publicprov.add_options()
        ("public-sharing-begins,G", value(s_.public_sharing_begins), "The period in which public sharing becomes available")
        ("public-sharing-tax,g", min<0>(s_.public_sharing_tax), "The per-period, lump sum tax collected from each reader for public sharing")
        ("prior-scale-public-sharing,S", min<1>(s_.prior_scale_public_sharing), "The same as --prior-weight, but applied in the first public sharing period")
        ;
    options_.add(publicprov);

    // This one is an object variable because the cli/gui need to add to it
    sim_controls_.add_options()
        ("periods,T", value(periods), "Number of simulation periods to run.")
        ("seed", value(seed), "Random seed to use.  If omitted, a random seed is obtained from the operating system's random source.")
        ("threads,j", value(threads), "Maximum number of threads to use for the simulation.  0 (the default) disables simulation threading entirely.")
        ;
    // Don't add this here: the caller has to do that (after adding to it, if necessary): boost
    // *copies* the argument to add(), so we can't add and then change it later.
    //options_.add(sim_options_);
}

CmdArgs::CLI::CLI(CreativitySettings &s) : Simulation(s) {}
CmdArgs::GUI::GUI(CreativitySettings &s) : Simulation(s) {}

// Single-letter options used:
// b B c C d D e f g G i j k K m M n N o O p P Q r R s S T u U w W x y z Z
// Available:
// a A E F H I J l L q t v V X Y

void CmdArgs::CLI::addOptions() {
    output = "creativity-SEED.crstate";
    threads = 0;

    Simulation::addOptions();

    sim_controls_.add_options()
        ("output,o", value(output), "Output file for simulation results.  If this contains the characters 'SEED', they will be replaced with the random seed value used for the simulation.")
        ("tmpdir", value(tmpdir), "Output directory in which to write the output file while running the simulation.  When "
            "the simulation finishes, the temporary file is moved to the output location specified by -o.  If this argument is omitted, the "
            "file is written directly to its final location.")
        ("overwrite,O", "Allows output file given to -o to be overwritten.")
        ;
    options_.add(sim_controls_);
}

void CmdArgs::GUI::addOptions() {
    threads = std::thread::hardware_concurrency();
    if (threads == 1) threads = 0;

    Simulation::addOptions();

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

void CmdArgs::Simulation::postParse(boost::program_options::variables_map &) {
    // If the user didn't give a seed, .seed won't have changed from seed, but we still don't want
    // to set it because explicitly setting a seed resets the RNG.
    if (eris::Random::seed() != seed) {
        eris::Random::seed(seed);
    }

    if (not output.empty()) {
        output = std::regex_replace(output, std::regex("SEED"), std::to_string(eris::Random::seed()));
    }

    s_.boundary = Creativity::boundaryFromDensity(s_.readers, s_.dimensions, density_);
}

void CmdArgs::CLI::postParse(boost::program_options::variables_map &vars) {
    Simulation::postParse(vars);
    overwrite = vars.count("overwrite") > 0;
}

void CmdArgs::GUI::postParse(boost::program_options::variables_map &vars) {
    Simulation::postParse(vars);
    start = vars.count("start") > 0;
    initialize = vars.count("initialize") > 0;
}


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
}
