#pragma once
#include "creativity/cmdargs/Simulator.hpp"
#include <string>

namespace boost { namespace program_options { class variables_map; } }
namespace creativity { struct CreativitySettings; }

namespace creativity { namespace cmdargs {

/** CmdArgs/Simulator subclass for graphical simulator command-line arguments. */
class GUI : public Simulator {
    public:
        /// Constructor for gui simulation arguments; takes the settings object
        GUI(CreativitySettings &s);

        /// Whether to start running the simulation in the GUI right away
        bool start = false;

        /** Whether to initialize (but not start) the simulation in the GUI right away.  Has no
         * effect if `start` is true.
         */
        bool initialize = false;

        /// The input file (for the GUI)
        std::string input;

    protected:
        /// Adds GUI command-line options into the option descriptions
        virtual void addOptions() override;

        /// Overridden to handle --initialize and --start options
        virtual void postParse(boost::program_options::variables_map &vars) override;

};

}}
