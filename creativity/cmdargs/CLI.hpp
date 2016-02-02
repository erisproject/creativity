#pragma once
#include "creativity/cmdargs/Simulator.hpp"
#include <string>

namespace boost { namespace program_options { class variables_map; } }
namespace creativity { struct CreativitySettings; }

namespace creativity { namespace cmdargs {

/** CmdArgs subclass for command-line simulator command-line arguments. */
class CLI : public Simulator {
    public:
        /// Constructor for cli simulation arguments; takes the settings object
        CLI(CreativitySettings &s);

        /// Set to tell cli to shut up (typically when running in batch mode)
        bool quiet = false;

        /// The temporary directory for results
        std::string tmpdir;

        /// Whether `output' can be overwritten
        bool overwrite = false;

        /// Overridden to add " -- command-line simulator"
        virtual std::string versionSuffix() const override;

    protected:
        /// Adds CLI command-line options into the option descriptions
        virtual void addOptions() override;

        /// Overridden to handle --overwrite option
        virtual void postParse(boost::program_options::variables_map &vars) override;

};

}}
