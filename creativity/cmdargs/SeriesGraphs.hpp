#pragma once
#include "creativity/cmdargs/CmdArgs.hpp"
#include "creativity/data/graph/Series.hpp"
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

        /// If true, we're in per-t confidence mode.  If false, we're in per-source file confidence mode.
        bool per_t_confidence = true;

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

        /** The title for the graph, with possible substitutions:
         *
         * %F - the full path of the input  as given in `input`
         * %p - the dirname of the input file (empty if no directory at all)
         * %f - the basename (i.e. filename) of the input path
         * %u - the unique portion of the input path
         * %w - the last word ([0-9A-Za-z_]) of the input file, not counting the extension.
         * %% - a literal %
         *
         * Any other %x string is reserved for future use.
         */
        std::string title = "%w (%u)";

        /// The width of the resulting graph, in inches.
        double width = 6;

        /// The height of the resulting graph, in inches.
        double height = 4;

        /** The `t` value at which to start the graph display.  If equal to -1, the graph display
         * range is determined by the input data.  Note that this only controls the graph range: the
         * observations can be further limited by setting data_min_t.
         */
        int graph_min_t = -1;

        /** The `t` value at which to end the graph display.  If equal to -1, the graph display
         * range is determined by the input data.  Note that this only controls the graph range:
         * observations can be further limited by setting data_max_t.
         */
        int graph_max_t = -1;

        /** The minimum t value required to actually include an observation in the graph.  -1 means
         * unlimited.
         *
         * The default is 1, since t=0 observations are pre-simulation observations and typically
         * irrelevant.
         * 
         * If this is less than graph_min_t it has no effect; if larger, observations between
         * graph_min_t and this won't be plotted (even though there is space for them on the graph).
         */
        int data_min_t = 1;

        /** The maximum t value to display in the graph.  -1 means unlimited.
         *
         * If this is greater than graph_max_t it has no effect; if less, observations between this
         * and graph_max_t won't be plotted (even though there is space for them on the graph).
         */
        int data_max_t = -1;

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

        /** The legend position */
        creativity::data::graph::Series::LegendPosition legend_position = creativity::data::graph::Series::LegendPosition::Right;

        /** The legend relative x value.
         * \sa creativity::data::graph::Series::legend_rel_x
         */
        double legend_rel_x = 1;

        /** The legend relative y value.
         * \sa creativity::data::graph::Series::legend_rel_y
         */
        double legend_rel_y = 0;

        /** The input string, such as "top-left", "outside-center", etc.  If blank, we just leave
         * legend_position/legend_rel_x/legend_rel_y alone. */
        std::string legend_position_input;

        /** The resolution (pixels per inch) for PNG output.  The PNG will have a pixel size of
         * this value times `width` and `height`, but will have the resolution encoded (so that
         * attempting to print it at default sizes should respect the width and height values).
         *
         * For PDF and SVG output, this value is ignored.
         */
        double resolution = 96;

        /// Overridden to handle t-min/t-max options
        virtual void postParse(boost::program_options::variables_map &vars) override;

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

