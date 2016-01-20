#pragma once
#include "creativity/data/graph/style.hpp"
#include "creativity/data/graph/Target.hpp"
#include <cairomm/context.h>
#include <map>

namespace creativity { namespace data { namespace graph {

/** Class to plot graphs of time series values and/or regions. */
class Series {
    public:
        /// No default constructor
        Series() = delete;

        /** Constructs a new series plot that plots to the given Cairo context.
         *
         * \param target the Target object to which the graph should be drawn
         * \param tmin the left-hand-side t value
         * \param tmax the right-hand-side t value
         * \param ymin the minimum y value at the bottom of the plot region
         * \param ymax the maximum y value at the top of the plot region
         */
        Series(Target &target, int tmin, int tmax, double ymin, double ymax);

        /** Writes the current graph to disk and begins a new graph.  If a multi-page surface
         * format, this will be a new page in the output file; otherwise it will be a new file.
         */
        void newPage();

        /** Plots the pairs of coordinates in the given map to the plot.  Keys should be the `t`
         * values.  If any sequential t values are missing or have non-finite values, the line will
         * be broken around the missing values.
         *
         * \param points the map of time to value coordinates to plot
         * \param style the LineStyle; defaults to a 1-unit-wide black line.
         */
        void addLine(const std::map<unsigned, double> &points,
                const LineStyle &style = LineStyle(RGBA(0,0,0))
                );

        /** Takes a map of time value to double pairs and plots the region between the pair values.
         * If any t values are missing or have one or more non-finite pair values, the region is
         * broken at that t value.
         *
         * \param intervals the map of time to region boundaries.  The order of the pair values does
         * not matter.
         * \param style FillStyle the fill style of the region.  The border of the fill style will
         * be used only for the top and bottom border of the region (i.e. the sides will have no
         * border).  Defaults to a 1/3 opacity blue fill with no (i.e. fully transparent) border.
         */
        void addRegion(const std::map<unsigned, std::pair<double, double>> &intervals,
                const FillStyle &style = FillStyle(RGBA(0,0,1,0.33333)));

    private:
        Target &target_;
        int tmin_, tmax_;
        double ymin_, ymax_;
        Cairo::Matrix translate_unit_, ///< Translates from [0,1]x[0,1] into surface space
            translate_graph_; ///< Translates from graph coordinates into surface space

        // The graph area boundary in the [0,1]x[0,1] graph size:
        static constexpr double graph_left_ = 0.1, graph_right_ = 0.8, graph_top_ = 0.1, graph_bottom_ = 0.8;

        // The size of the graph border in surface units.
        double graph_border_ = 1;

        // The buffer to add on the top/bottom edge of the graph (so that lines/points don't run
        // onto the graph outline).  In absolute surface units.
        double graph_buffer_tb_ = 3;
        // The buffer to add to the left/right edges.  Defaults to 0 (because L/R values won't run
        // over anyway).
        double graph_buffer_lr_ = 0;

        // Whether initializePage() has been called for the current page yet
        bool page_initialized_ = false;

    protected:
        /** Draws the initial layout (graph box, etc.).  Does nothing if called repeatedly (until
         * the next newPage() call).
         */
        void initializePage();

        /** Clips the given Cairo::Context to the graph space.  This should be called without a
         * transformation having been applied to the surface.
         */
        void clipToGraph(Cairo::RefPtr<Cairo::Context> ctx) const;

        /** Helper method called by addRegion: takes an iterator to the first and final elements of
         * a contiguous set, draws a region empassing the min/max values, and closes it so that a
         * subsequent fill will fill it in.  If the region consists of only a single point (that is,
         * first == last), the point is copied into strays to be later drawn as vertical lines by
         * addRegion().
         *
         * \param ctx the Cairo::Context, which should have a transformation applied so that data
         * points correspond to the graph region.
         * \param first the iterator to the first element in the group
         * \param last the iterator to the last element in the group (NOT the past-the-end iterator)
         * \param stays a map to store a "stray" (that is, single-time-period) observation in (only
         * when first == last).
         */
        static void boundRegion(
                Cairo::RefPtr<Cairo::Context> ctx,
                std::map<unsigned, std::pair<double, double>>::const_iterator first,
                std::map<unsigned, std::pair<double, double>>::const_iterator last,
                std::map<unsigned, std::pair<double, double>> &strays);
};

}}}
