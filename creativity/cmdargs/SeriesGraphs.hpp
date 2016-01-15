#pragma once
#include "creativity/cmdargs/CmdArgs.hpp"
#include <limits>
#include <string>
#include <vector>
#include <set>

namespace creativity { namespace cmdargs {

/** CmdArgs subclass for creativity-series-plots arguments (for plotting quantiles from series
 * and/or quantiles files). */
class SeriesGraphs : public CmdArgs {
    public:
        /// Constructor for series arguments; takes no arguments
        SeriesGraphs();

        /** The confidence levels to plot as a comma-separated list.  If input is a series file, the
         * needed quantiles are calculated (as would be done by creativity-series-quantiles); if given a
         * quantiles file, the necessary quantiles must exist in the file.
         *
         * A value of .9, for example, uses the .05 and .95 quantiles.
         *
         * The special value of 0 corresponds to the median.
         *
         * Order is not important, and duplicate values are ignored.
         *
         * The default plots the median, 50%, 90%, and 95% confidence intervals.
         */
        std::string levels = "median,0.5,0.9,0.95";

        /** The input files (series or quantiles) to plot.  Plot axes scales will be identical
         * across files, so typically this should be called with series files for the same (or at
         * least directly comparable) variables.
         */
        std::vector<std::string> input;

        /** The output filename, which must end in ".pdf", ".png", or ".svg".  If this ends in
         * ".pdf" or ".svg", the file will be a single PDF or SVG file with one page per input file;
         * if the output file is "something.png", files "something-1.png", ..., "something-n.png",
         * one image per input file.
         *
         * The default is "series-plots.pdf".
         */
        std::string output = "series-plots.pdf";

        /** Whether to overwrite the output file(s) if it (they) exist.
         */
        bool overwrite = false;

        /// The width of the resulting graph, in inches.
        double width = 6;

        /// The height of the resulting graph, in inches.
        double height = 4;

        /** The resolution (pixels per inch) for PNG output.  The PNG will have a pixel size of
         * this value times `width` and `height`, but will have the resolution encoded (so that
         * attempting to print it at default sizes should respect the width and height values).
         *
         * For PDF and SVG output, this value is ignored.
         */
        double resolution = 96;

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

