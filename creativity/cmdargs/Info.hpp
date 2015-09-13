#pragma once
#include "creativity/cmdargs/CmdArgs.hpp"
#include <eris/types.hpp>
#include <limits>
#include <string>
#include <vector>

namespace boost { namespace program_options { class variables_map; } }

namespace creativity { namespace cmdargs {

/** CmdArgs subclass for creativity-info script to display simulation information. */
class Info : public CmdArgs {
    public:
        /// The number of significant digits in output values
        unsigned int output_precision = 6;

        /// The input file to load data from
        std::string input;

        /// Overridden to add " FILENAME"
        virtual std::string usage() const override;

        /// Overridden to add info about the input file
        virtual std::string help() const override;

        /// Overridden to add " -- simulation detail viewer"
        virtual std::string versionSuffix() const override;

    protected:
        /// Adds data collector command-line options into the option descriptions
        virtual void addOptions() override;

        /// Overridden to make sure the input file is specified and non-empty
        virtual void postParse(boost::program_options::variables_map &vars) override;

};

}}
