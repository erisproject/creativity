#pragma once
#include "creativity/cmdargs/CmdArgs.hpp"
#include <eris/types.hpp>
#include <limits>
#include <string>
#include <vector>

namespace boost { namespace program_options { class variables_map; } }

namespace creativity { namespace cmdargs {

/** CmdArgs/Simulator subclass for graphical simulator command-line arguments. */
class Data : public CmdArgs {
    public:
        /// Constructor for gui simulation arguments; takes the settings object
        Data();

        /// Struct for storing --verify-* options
        struct {
            /// If non-zero, this verifies that given crstate data files have the same piracy_begins parameter
            eris::eris_time_t piracy_begins = 0;
            /// If non-zero, verify the total number of simulation periods
            eris::eris_time_t periods = 0;
            /// If non-zero, verify the public_sharing_begins of given simulations
            eris::eris_time_t public_sharing_begins = 0;
        } verify;

        /// Struct for storing --skip-* options
        struct {
            /// Skip piracy data if true
            bool piracy = false;
            /// Skip public sharing data if true
            bool public_sharing = false;
        } skip;

        /// If true, produce output in human-readable format; otherwise CSV output.
        bool human_readable = false;

        /// The number of periods to use to calculate average values
        unsigned int data_periods = 25;

        /// The precision of double values.  The default is full double precision.
        unsigned int double_precision = std::numeric_limits<double>::max_digits10;

        /// The input files to load data from
        std::vector<std::string> input;

        /// Overridden to add " CRSTATE [CRSTATE ...]"
        virtual std::string usage() const override;

        /// Overridden to add info about input files
        virtual std::string help() const override;

    protected:
        /// Adds data collector command-line options into the option descriptions
        virtual void addOptions() override;

        /// Overridden to handle --initialize and --start options
        virtual void postParse(boost::program_options::variables_map &vars) override;

};

}}
