#include "creativity/cmdargs/Simulator.hpp"
#include "creativity/cmdargs/strings.hpp"
#include <eris/types.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/errors.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <cstdint>
#include <regex>
#include <sstream>
namespace boost { namespace program_options { class variables_map; } }


namespace creativity { namespace cmdargs {

namespace po = boost::program_options;

Simulator::Simulator(CreativitySettings &cs) : s_(cs) {}

void Simulator::addOptions() {
    CmdArgs::addOptions();
    po::options_description structure("Structure"), initial("Initial Behaviour"),
        authorship("Authorship Settings and Costs"), costs("Costs"), beliefs("Beliefs"),
        piracy("Piracy"), publicprov("Public Provisioning");

    structure.add_options()
        ("dimensions,D", min<1>(s_.dimensions), "    Number of dimensions of the simulation")
        ("readers,r", min<1>(s_.readers), "    Number of reader/author agents in the simulation")
        ("density,d", above<0>(density_), (u8"  Reader density (in readers per unit^D, where D is the configured # of dimensions).  The default is the density required to have simulation boundaries at ±" + output_string(s_.boundary) + " in each dimension.").c_str())
        ("reader-step-sd,R", min<0>(s_.reader_step_sd), "Standard deviation of the inter-period random-direction reader movement distance ~ |N(0, R)|")
        ("book-distance-sd,B", min<0>(s_.book_distance_sd), "Standard deviation of book distance from author; distance ~ |N(0, B)|")
        ("book-quality-sd,Q", min<0>(s_.book_quality_sd), "    Standard deviation of book perceived quality; perceived quality ~ N(q, Q), where q is the innate quality and Q is this value")
        ;
    options_.add(structure);

    authorship.add_options()
        ("creation-time,e", value(s_.creation_time), u8"The number of periods that elapse between the creation decision and the book being ready")
        ("reader-creation-shape,s", below<1>(s_.reader_creation_shape), u8"Shape parameter, β, of the creator effort function: q(ℓ) = (α/β)[(ℓ+1)^β - 1]")
        ("reader-creation-scale-min,z", min<0>(s_.reader_creation_scale_min), u8"Minimum support of scale parameter α ~ U[a,b] for the creator effort function: q(ℓ) = (α/β)[(ℓ+1)^β - 1]")
        ("reader-creation-scale-max,Z", min<0>(s_.reader_creation_scale_max), u8"Maximum support of scale parameter α ~ U[a,b] for the creator effort function: q(ℓ) = (α/β)[(ℓ+1)^β - 1]")
        ;
    options_.add(authorship);

    costs.add_options()
        ("income,i", above<0>(s_.income), "  Per-period external reader income")
        ("cost-fixed,C", min<0>(s_.cost_fixed), "    Fixed cost of keeping a book on the market for a period")
        ("cost-unit,c", min<0>(s_.cost_unit), "    Unit cost of making a copy of a book")
        ("cost-piracy,y", min<0>(s_.cost_piracy), "    Cost of receiving a pirated copy of a book")
        ;
    options_.add(costs);

    beliefs.add_options()
        ("prior-scale,w", min<1>(s_.prior_scale), "    The per-period standard deviation scaling factor applied when a previous belief is used as the next period's prior")
        ("prediction-draws,p", min<1>(s_.prediction_draws), "The number of draws agents take from beliefs when using for prediction.  Lower values can be used to introduce deliberate randomness into the simulation")
        ("burnin-periods,u", value(s_.burnin_periods), "    The number of initial periods during which `--prior-scale-burnin' should be used instead of `--prior-scale'")
        ("belief-threshold,b", value(s_.initial.belief_threshold), "  The minimum n-k value at which a readers bases decision on beliefs instead of initial parameters")
        ("prior-scale-burnin,U", min<1>(s_.prior_scale_burnin), "The same as --prior-weight, but applied in the first `--burnin-periods' periods")
        ;
    options_.add(beliefs);

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

    piracy.add_options()
        ("piracy-begins,P", value(s_.piracy_begins), "    The period in which piracy becomes available.  0 means never")
        ("piracy-link-proportion,f", range<0, 1>(s_.piracy_link_proportion), "Proportion of potential sharing links between readers that are created")
        ("prior-scale-piracy,W", min<1>(s_.prior_scale_piracy), "The same as --prior-weight, but applied in the first piracy period")
        ;
    options_.add(piracy);

    publicprov.add_options()
        ("public-sharing-begins,G", value(s_.public_sharing_begins), "The period in which public sharing becomes available. 0 means never")
        ("public-sharing-tax,g", min<0>(s_.public_sharing_tax), "The per-period, lump sum tax collected from each reader for public sharing")
        ("prior-scale-public-sharing,S", min<1>(s_.prior_scale_public_sharing), "The same as --prior-weight, but applied in the first public sharing period")
        ;
    options_.add(publicprov);

    // This one is an object variable because the cli/gui need to add to it
    sim_controls_.add_options()
        ("periods,T", value(periods), "    Number of simulation periods to run.")
        ("seed", value(seed), "    Random seed to use.  If omitted, a random seed is obtained from the operating system's random source.")
        ("threads,j", value(threads), "    Maximum number of threads to use for the simulation.  0 (the default) disables simulation threading entirely.")
        ;
    // Don't add this here: the caller has to do that (after adding to it, if necessary): boost
    // *copies* the argument to add(), so we can't add and then change it later.
    //options_.add(sim_options_);
}


void Simulator::postParse(boost::program_options::variables_map &) {
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

}}
