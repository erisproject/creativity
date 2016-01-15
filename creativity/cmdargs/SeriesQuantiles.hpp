#pragma once
#include "creativity/cmdargs/CmdArgs.hpp"
#include <limits>
#include <string>
#include <vector>
#include <set>

namespace creativity { namespace cmdargs {

/** CmdArgs subclass for creativity-series-quantiles arguments (for extracting quantiles from series files). */
class SeriesQuantiles : public CmdArgs {
    public:
        /// Constructor for series arguments; takes no arguments
        SeriesQuantiles();

        /** Quantiles to extract.  If not specified, defaults to {0, 0.005, 0.01, 0.025, 0.05, 0.25,
         * 0.5, 0.75, 0.95, 0.975, 0.99, 0.995, 1}.  The 0 quantile is defined to be the minimum
         * value, the 1 quantile the maximum.
         *
         * When a requested quantile lies between observations, linear interpolation between
         * immediate neighbours is used to construct the quantile.  For example, if the data
         * contains 5 values, {0, 1, 2, 5, 100}, the 0, 0.25, 0.5, 0.75, and 1 quantiles match the
         * values 0, 1, 2, 5, and 100 exactly.  The 0.95 quantile, for example, will equal 0.2(5) +
         * 0.8(100) = 81.
         *
         * The resulting file will have header values qXXX where XXX are the the significant digits
         * following the decimal place, e.g. q005, q75, q975 for .005, .75, and .975 quantiles.  The
         * special values 0, 0.5, and 1 are replaced with min, median, and max, respectively.
         */
        std::string quantiles = "0,.005,.01,.025,.05,.25,.5,.75,.95,.975,.99,.995,1";

        /// The precision of double values.  The default is 10.
        unsigned int double_precision = 10;

        /// The input files to load data from
        std::vector<std::string> input;

        /** Remove this filename prefix, if found.  If an input file starts with the given string,
         * it is removed when constructing the output filename.
         *
         * Defaults to "series-".
         *
         * \sa output_prefix
         */
        std::string output_unprefix = "series-";

        /** The output prefix; output files will be placed in the same directory as input files,
         * with this string prepended to the beginning of the filename.  The prefix is added
         * removing the `output_strip` prefix (if found).
         *
         * Defaults to "quantiles-".
         *
         * \sa output_strip
         */
        std::string output_prefix = "quantiles-";

        /** Whether to overwrite output files if they exist.  (Note that input files will never be
         * overwritten, even if this is true.)
         */
        bool overwrite = false;

        /// Overridden to add " FILE [FILE ...]"
        virtual std::string usage() const override;

        /// Overridden to add info about input files
        virtual std::string help() const override;

        /// Overridden to add " -- simulation data series generator"
        virtual std::string versionSuffix() const override;

    protected:
        /// Adds series command-line options into the option descriptions
        virtual void addOptions() override;
};

}}
