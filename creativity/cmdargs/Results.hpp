#pragma once
#include "creativity/cmdargs/CmdArgs.hpp"
#include <eris/types.hpp>
#include <limits>
#include <string>
#include <vector>

namespace boost { namespace program_options { class variables_map; } }

namespace creativity { namespace cmdargs {

/** CmdArgs subclass for creativity-results script to run models on generated data. */
class Results : public CmdArgs {
    public:
        /// The supported format types
        enum class Format {
            Text, ///< enum value for plain-text output
            HTML, ///< enum value for HTML output
            LaTeX ///< enum value for LaTeX output
        };

        /// The requested format type
        Format format = Format::Text;

        /// Only one of these should be set to set the format
        bool format_text = false, format_html = false, format_latex = false;

        /// The number of significant digits in output values
        unsigned int output_precision = 6;

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
