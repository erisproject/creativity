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

        /// The `t` value at which to start the graph display
        unsigned graph_min_t = 0;

        /** The `t` value at which to end the graph display.  If 0, the value graph display end is
         * determined by the data_stretch value.
         */
        unsigned graph_max_t = 0;

        /** The `t` value at which to start plotting data; has no effect if less than graphmin, or
         * less than the first observation with data. */
        unsigned data_min_t = 0;

        /** Determines the minimum value to include in the series plot.  The default, NaN,
         * determines the minimum from the data.
         */
        double graph_min_value = std::numeric_limits<double>::quiet_NaN();

        /** Determines the maximum value to include in the series plot.  The default, NaN,
         * determines the maximum from the data.
         */
        double graph_max_value = std::numeric_limits<double>::quiet_NaN();

        /** Determines what happens when input series have different maximum time periods, and
         * graph_max_t is set to 0.  If false, the horizontal axes of graphs will be scaled so that
         * all graphs fill the graph space, but axes are not directly comparable--that is, each
         * graph has a maximum graph `t` value determined by its file.  If true, the
         * implicit graph maximum will be determined by the maximum `t` value across all input
         * files, not just the one being graphed.
         *
         * For example, if called with two input files A.csv and B.csv, where A.csv has periods t=0
         * through t=100 and B.csv has t=0 through t=200, when this is false, the graph for A will
         * go from 0 to 100 while B's goes from 0 to 200; when this is true, both graphs will run
         * from 0 to 200 (and the plotted data in A will simply stop at 100).
         *
         * The default is false (adjust the axes so each plot takes up the full width).
         */
        bool same_horizontal_scale = false;

        /** Controls vertical axes when graph_min_value or graph_max_value are not set to finite
         * values.  When false the vertical axis of each plot will be set so as to include all
         * plotted quantile values for only the values in the displayed file.  If true, the vertical
         * axes will be set to include all displayed quantile values for all displayed files, thus
         * making the vertical axis directly comparable across plots.
         *
         * The default is true (axes are the same across graphs).
         */
        bool same_vertical_scale = true;

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

