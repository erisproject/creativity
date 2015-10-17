#pragma once
#include "creativity/cmdargs/CmdArgs.hpp"
#include "creativity/data/tabulate.hpp"
#include <eris/types.hpp>
#include <limits>
#include <string>
#include <vector>

namespace boost { namespace program_options { class variables_map; } }

namespace creativity { namespace cmdargs {

/** CmdArgs subclass for creativity-results script to run models on generated data. */
class Results : public CmdArgs {
    public:
        /// The requested format data
        struct {
            /// The output format type, defaults to text
            data::TableFormat type = data::TableFormat::Text;

            /// Whether --text was explicitly requested
            bool text = false;
            /// Whether --html was explicitly requested
            bool html = false;
            /// Whether --latex was explicitly requested
            bool latex = false;

            /// The number of significant digits in output values
            unsigned int precision = 6;
        } format;

        /// Whether to run various types of analysis
        struct {
            /// Conditional analysis of parameter distributions in write vs non-write simulations
            bool write_or_not = false;
            /// Include correlation/covariance matrix in write_or_not analysis
            bool write_or_not_corrcov = false;
            /// Average effects model
            bool average = false;
            /// Marginal effects model
            bool marginal = false;
            /** All of the above; this does not need to be checked: if specified, all of the above
             * will be set to true except for write_or_not_corrcov, which needs to be specified
             * explicitly.
             */
            bool all = false;
            /// None of the above; none of the above will be true.
            bool none = false;
        } analysis;

        /// Output files
        struct {
            std::string filename; ///< Result analysis output; empty means output to stdout

            std::string dump_writing; ///< CSV output of "writing" data
            std::string dump_nowriting; ///< CSV output of "nowriting" data

            bool overwrite = false; ///< If true, overwrite the above, otherwise throw exception.
        } output;


        /// The input file to load data from
        std::string input;

        /// Overridden to add " FILENAME"
        virtual std::string usage() const override;

        /// Overridden to add info about the input file
        virtual std::string help() const override;

        /// Overridden to add " -- result analyzer"
        virtual std::string versionSuffix() const override;

    protected:
        /// Adds data collector command-line options into the option descriptions
        virtual void addOptions() override;

        /// Overridden to make sure only one of --text/--html/--latex are given.
        virtual void postParse(boost::program_options::variables_map &vars) override;

};

}}
