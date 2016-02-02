#pragma once
#include "creativity/cmdargs/CmdArgs.hpp"
#include "creativity/Creativity.hpp"
#include <boost/program_options/options_description.hpp>
#include <eris/Random.hpp>
#include <string>
namespace boost { namespace program_options { class variables_map; } }

namespace creativity { namespace cmdargs {

/** Common base class for CLI and GUI argument handling for running a simulation.
 *
 * Current used single-argument parameters:
 *
 * Single-letter options used in this class:
 *     b B c C d D e E f g G i j k K m M n N p P Q r R s S T u U w W x y z Z
 * Used in [CLI], {GUI}, or both:
 *     o [O] [q]
 * Used in CmdArgs base class:
 *     h
 * Available:
 *     a A F H I J l L t v V X Y
 */
class Simulator : public CmdArgs {
    protected:
        /// Default constructor deleted
        Simulator() = delete;

        /// Constructs a Simulator object that stores values in the given CreativitySettings object.
        Simulator(CreativitySettings &cs);

    public:
        // Variables which can't be stored directly in CreativitySettings:
        /** The number of periods to run the simulation.  Before calling addCliOptions() or
         * addGuiOptions(), this is the default value; after calling parse this will be updated
         * to whatever the user specified.
         */
        unsigned int periods = 300;

        /** The output file for simulation results; default is "creativity-SEED.crstate" for the
         * CLI, unset for the GUI.
         */
        std::string output;

        /** The number of threads to use.  The default is number of hardware threads for the
         * GUI, 0 for the CLI.
         */
        unsigned int threads = 0;

        /** The seed.  The default is whatever eris::Random::seed() returns, which is random
         * (unless overridden with ERIS_RNG_SEED).  This value can be ignored: it is handled by
         * parse().
         */
        typename eris::Random::rng_t::result_type seed = eris::Random::seed();

    protected:
        /** Adds common options into the options descriptions.  Called by
         * CLI::addOptions() and GUI::addOptions().
         */
        void addOptions() override;

        /// Overridden to handle output, seed, and and boundary settings
        virtual void postParse(boost::program_options::variables_map &vars) override;

        /** The options (which will be added to the end of options_) containing simulator settings,
         * some of which differ between the CLI and the GUI.
         */
        boost::program_options::options_description sim_controls_{"Simulator controls"};

        /** The settings reference in which to store given values */
        CreativitySettings &s_;

        /// Density pseudo-parameter: the stored value is actually boundary
        double density_ = Creativity::densityFromBoundary(s_.readers, s_.dimensions, s_.boundary);

};

}}
