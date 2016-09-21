#include "creativity/cmdargs/Simulator.hpp"
#include "creativity/cmdargs/strings.hpp"
#include <eris/types.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/errors.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/variables_map.hpp>
#include <cstdint>
#include <regex>
#include <sstream>


namespace creativity { namespace cmdargs {

namespace po = boost::program_options;

Simulator::Simulator(CreativitySettings &cs) : s_(cs) {}

void Simulator::addOptions() {
    CmdArgs::addOptions();
    po::options_description structure("Structure"), initial("Initial Behaviour"),
        authorship("Authorship Settings and Costs"), costs("Costs"), beliefs("Beliefs"),
        piracy("Piracy"), policy("Policy options"), pol_public("Public sharing policy (requires --policy=public-sharing)"),
        pol_catch("Catch and fine policy (requires --policy=catch-pirates)");

    structure.add_options()
        ("dimensions,D", min<1>(s_.dimensions), "    Number of dimensions of the simulation")
        ("readers,r", min<1>(s_.readers), "    Number of reader/author agents in the simulation")
        ("density,d", above<0>(density_), (u8"  Reader density (in readers per unit^D, where D is the configured # of dimensions).  The default is the density required to have simulation boundaries at ±" + output_string(s_.boundary) + " in each dimension.").c_str())
        ("reader-step-mean,R", min<0>(s_.reader_step_mean), "Mean of the inter-period random-direction reader movement distance.  The distance stepped is this value times a draw from a χ² with 1 d.f.")
        ("book-distance-mean,B", min<0>(s_.book_distance_mean), "Mean of book distance from author.  The value is this times a draw from a χ² with 1 d.f.")
        ("book-quality-sd,Q", min<0>(s_.book_quality_sd), "    Standard deviation of book perceived quality; perceived quality ~ N(q, Q), where q is the innate quality and Q is this value")
        ;
    options_.add(structure);

    authorship.add_options()
        ("creation-fixed,E", value(s_.creation_fixed), u8"The fixed cost component of creation; total initial creation cost is this plus expended effort")
        ("creation-time,e", value(s_.creation_time), u8"The number of periods that elapse between the creation decision and the book being ready")
        ("reader-creation-shape,s", below<1>(s_.reader_creation_shape), u8"Shape parameter, β, of the creator effort function: q(ℓ) = (α/β)[(ℓ+1)^β - 1]")
        ("reader-creation-scale-min,z", min<0>(s_.reader_creation_scale_min), u8"Minimum support of scale parameter α ~ U[a,b] for the creator effort function: q(ℓ) = (α/β)[(ℓ+1)^β - 1]")
        ("reader-creation-scale-range,Z", min<0>(s_.reader_creation_scale_range), u8"Support range (that is, b-a) of scale parameter α ~ U[a,b] for the creator effort function: q(ℓ) = (α/β)[(ℓ+1)^β - 1]")
        ;
    options_.add(authorship);

    costs.add_options()
        ("income,i", above<0>(s_.income), "  Per-period external reader income")
        ("cost-market,C", min<0>(s_.cost_market), "    Fixed cost of keeping a book on the market for a period")
        ("cost-unit,c", min<0>(s_.cost_unit), "    Unit cost of making a copy of a book")
        ("cost-piracy,y", min<0>(s_.cost_piracy), "    Cost of receiving a pirated copy of a book")
        ;
    options_.add(costs);

    beliefs.add_options()
        ("prior-scale,w", min<1>(s_.prior_scale), "    The per-period standard deviation scaling factor applied when a previous belief is used as the next period's prior")
        ("prediction-draws,p", min<1>(s_.prediction_draws), "The number of draws agents take from beliefs when using for prediction.  Lower values can be used to introduce deliberate randomness into the simulation")
        ("burnin-periods,u", value(s_.burnin_periods), "    The number of initial periods during which `--prior-scale-burnin' should be used instead of `--prior-scale'")
        ("belief-threshold,b", value(s_.initial.belief_threshold), "  The minimum n-k value at which a readers bases decision on beliefs instead of initial parameters")
        ("prior-scale-burnin,U", min<1>(s_.prior_scale_burnin), "The same as --prior-scale, but applied in the first `--burnin-periods' periods")
        ;
    options_.add(beliefs);

    initial.add_options()
        ("initial-prob-write,J", range<0, 1>(s_.initial.prob_write), "The probability of writing in initial periods")
        ("initial-effort-min,m", min<0>(s_.initial.l_min), "The minimum support of effort l ~ U[a,b] for authored books in initial periods")
        ("initial-effort-range,M", min<0>(s_.initial.l_range), "The range of effort l ~ U[a,b] (in other words, b-a) for authored books in initial periods")
        ("initial-price-min,n", min<0>(s_.initial.p_min), "The minimum support, `a', of the random component of price c + U[a,b] for new books in initial periods")
        ("initial-price-range,N", min<0>(s_.initial.p_range), "The range, `b-a', of the random component of price c + U[a,b] for new books in initial periods")
        ("initial-prob-keep,k", range<0, 1>(s_.initial.prob_keep), "The probability of keeping a previously-written book on the market for another period")
        ("initial-keep-price,K", range<0, 1>(s_.initial.keep_price), "The  price-above-marginal-cost level (relative to current P-MC) for a book left on the market for another period")
        ;
    options_.add(initial);

    piracy.add_options()
        ("piracy-begins,P", value(s_.piracy_begins), "    The period in which piracy becomes available.  0 means never")
        ("piracy-link-proportion,L", range<0, 1>(s_.piracy_link_proportion), "Proportion of potential sharing links between readers that are created")
        ("prior-scale-piracy,W", min<1>(s_.prior_scale_piracy), "The same as --prior-scale, but applied in the first piracy period")
        ;
    options_.add(piracy);

    policy.add_options()
        ("policy,g", value(policies), "A list of policies to enable in the policy phase of the simulation.  This option may be specified multiple times, or may be specified with a comma-separated list of policies.  If the option is omitted, defaults to `public-sharing`; specify 'none' to disable all policies.  Currently accepted policy names: 'public-sharing' - enables public sharing with redistribution (can be shorted to 'public'); 'catch-pirates' - enables detection and fines (can be shorted to 'catch')")
        ("policy-begins,G", value(s_.policy_begins), "When a policy response is enabled, this specifies the period in which it begins.  Has no effect if no policy response is enabled")
        ("prior-scale-policy,S", min<1>(s_.prior_scale_policy), "The same as --prior-scale, but applied in the first policy response period")
        ;

    pol_public.add_options()
        ("public-sharing-tax,A", min<0>(s_.policy_public_sharing_tax), "The per-period, lump sum tax collected from each reader for public sharing")
        ;
    policy.add(pol_public);

    pol_catch.add_options()
        ("catch-tax,x", min<0>(s_.policy_catch_tax), "The lump-sum, per-reader tax that is collected to pay for the catch policy")
        ("catch-cost", min<0>(s_.policy_catch_cost), "The cost that is incurred by someone accused of piracy, even if innocent")
        ("catch-fine-const", value(s_.policy_catch_fine[0]), "The constant term, c, of the fine polynomial, c + bP + aP²")
        ("catch-fine-lin", value(s_.policy_catch_fine[1]), "The linear term coefficient, b, of the fine polynomial, c + bP + aP²")
        ("catch-fine-sq", value(s_.policy_catch_fine[2]), "The squared term coefficient, a, of the fine polynomial, c + bP + aP²")
        ("catch-mu-const", value(s_.policy_catch_mu[0]), "The constant term, c, that determines μ = c + bF + aF², where F is the --catch-tax value")
        ("catch-mu-lin", value(s_.policy_catch_mu[1]), "The linear term coefficient, b, that determines μ = c + bF + aF², where F is the --catch-tax value")
        ("catch-mu-sq", value(s_.policy_catch_mu[2]), "The squared term coefficient, a, that determines μ = c + bF + aF², where F is the --catch-tax value")
        ("catch-sigma-const", value(s_.policy_catch_sigma[0]), "The constant term, c, that determines σ = c + bF + aF², where F is the --catch-tax value")
        ("catch-sigma-lin", value(s_.policy_catch_sigma[1]), "The linear term coefficient, b, that determines σ = c + bF + aF², where F is the --catch-tax value")
        ("catch-sigma-sq", value(s_.policy_catch_sigma[2]), "The squared term coefficient, a, that determines σ = c + bF + aF², where F is the --catch-tax value")
        ;
    policy.add(pol_catch);

    options_.add(policy);

    // This one is an object variable because the cli/gui need to add to it
    sim_controls_.add_options()
        ("periods,T", value(periods), "    Number of simulation periods to run.")
        ("seed", value(seed), "    Random seed to use.  If omitted, a random seed is obtained from the operating system's random source.")
        ("threads,j", value(threads), "    Maximum number of threads to use for the simulation.  0 (the default) disables simulation threading entirely.")
        ("tmpdir", value(tmpdir), "Output directory in which to write the output file while running the simulation.  When "
            "the simulation finishes, the temporary file is compressed and the compression version written to the final destination "
            "specified by -o.  If this argument is omitted, the file is written to a temporary file in the same directory as the final "
            "output file.  Has no effect when --memory is enabled.")
        ;
    // Don't add this here: the caller has to do that (after adding to it, if necessary): boost
    // *copies* the argument to add(), so we can't add and then change it later.
    //options_.add(sim_options_);
}


void Simulator::postParse(boost::program_options::variables_map &vars) {
    // If the user didn't give a seed, .seed won't have changed from seed, but we still don't want
    // to set it because explicitly setting a seed resets the RNG.
    if (eris::random::seed() != seed) {
        eris::random::seed(seed);
    }

    if (not output.empty()) {
        output = std::regex_replace(output, std::regex("SEED"), std::to_string(eris::random::seed()));
    }

    s_.boundary = Creativity::boundaryFromDensity(s_.readers, s_.dimensions, density_);

    // Default, but only if the argument isn't specified (so that `--policy=` by itself disables the
    // policy)
    if (vars.count("policy") == 0) {
        policies.push_back("public-sharing");
    }

    for (const auto &p : policies) {
        std::istringstream iss;
        iss.str(p);
        std::string policy;
        while (std::getline(iss, policy, ',')) {
            if (policy == "public" or policy == "public-sharing")
                s_.policy |= POLICY_PUBLIC_SHARING;
            else if (policy == "catch" or policy == "catch-pirates")
                s_.policy |= POLICY_CATCH_PIRATES;
            else if (policy == "none" or policy == "")
                /* ignore */;
            else
                throw po::invalid_option_value("--policy " + policy);
        }
    }
}

}}
