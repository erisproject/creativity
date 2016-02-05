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

        /// If true, decompress xz files into memory
        bool memory_xz = false;

        /// If true, copy all files into memory
        bool memory = false;

        /// Show every `n`th period (where this value is n).  0 means don't show any.
        unsigned int thin_periods = 10;

        /// Show CLI usage to recreate experiment
        bool show_cli_args = false;

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
