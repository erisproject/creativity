#pragma once
#include "creativity/data/graph/style.hpp"
#include "creativity/data/graph/Target.hpp"
#include <cairomm/context.h>
#include <pangomm/fontdescription.h>
#include <pangomm/font.h>
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
         * \param title_markup the title of the graph (which may contain pango markup)
         * \param tmin the left-hand-side t value
         * \param tmax the right-hand-side t value
         * \param ymin the minimum y value at the bottom of the plot region
         * \param ymax the maximum y value at the top of the plot region
         */
        Series(Target &target, std::string title_markup, int tmin, int tmax, double ymin, double ymax);

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
                const FillStyle &style = FillStyle(RGBA(0,0,1,1./3.)));

        /** Adds a simple legend item.  The legend box has two main drawing components; each can be
         * transparent to suppress the element.
         *
         * \param markup the text (with pango markup) to display for the legend item
         * \param lower fill style for the bottom half of the box.  The border component of the fill
         * style is drawn across the middle of the box.
         * \param preserve whether this legend item should be automatically added to new pages when
         * newPage() is called.
         * \param background the fill and border style for the legend box.  Defaults to a
         * transparent background with a 1-unit-wide black border.  Note that the background fill is
         * always contained entirely within the border, even if the border is transparent.
         *
         * The background fill and border is drawn first, followed by the lower area.
         */
        void addLegendItem(std::string markup,
                FillStyle lower,
                bool preserve = true,
                FillStyle background = FillStyle(Transparent, LineStyle(RGBA(0,0,0), 1)));

        /// Alias for the function callback for custom legend support
        using LegendPainterCallback_t = std::function<void(Cairo::RefPtr<Cairo::Context>, const Cairo::Matrix&)>;

        /** Adds a legend item with custom box drawing.  The custom callback function is called with
         * a Cairo::Context object to allow the caller to draw the box as desired.
         *
         * The `legend_painter` callback is called with the Cairo::Context object with an applied
         * clip region that allows drawing only in the legend box and legend text spaces.  Two
         * matrices are also supplied: one that translated [0,1]x[0,1] coordinates into the region
         * where the legend box should be drawn, and the other that translates [0,1]x[0,1]
         * coordinates into the legend text rectangular area.  Note that (0,0) and (1,1) correspond
         * to the upper-left and lower-right corners of the respective regions (so larger `y` values
         * are down the page, in the usual graphical toolkit drawing convention).
         *
         * \param name the legend name
         * \param legend_painter a callable function object that receives a Cairo::Context to draw to
         * and a Cairo::Matrix object mapping [0,1]x[0,1] coordinates into the (clipped) legend box.
         * \param preserve if true (the default), this legend item will be automatically re-added to
         * the graph whenever newPage() is called to start a new graph.  If false, the legend item
         * is only added to the current page.
         */
        void addLegendItem(std::string name, const LegendPainterCallback_t &legend_painter, bool preserve = true);

        /// Same as above, but takes an rvalue callback object.
        void addLegendItem(std::string name, LegendPainterCallback_t &&legend_painter, bool preserve = true);

        /** Draws a rectangle border, fill the rectangle interior, and (optionally) clip to the
         * rectangle interior.
         *
         * \param ctx a Cairo::Context positioned at the top-left corner of the desired rectangle
         * \param width the desired width (including the border) of the rectangle
         * \param height the desired height (including the border) of the rectangle
         * \param style a FillStyle from which to read the border style and background fill colour
         * \param clip if true, clip the context to the interior of the rectangle before returning
         */
        static void drawRectangle(Cairo::RefPtr<Cairo::Context> ctx, double width, double height, const FillStyle &style, bool clip = false);

        /** The surface units between the top graph border and top surface (image) edge. 
         * Note that this space contains the graph title.
         */
        double graph_top = 40;

        /// The surface units between the left graph border and left surface (image) edge
        double graph_left = 40;

        /// The surface units between the bottom graph border and bottom surface (image) edge
        double graph_bottom = 40;

        /** The surface units between the right graph border and right surface (image) edge.  Note
         * that this space typically contains any legend items.
         */
        double graph_right = 80;

        /** The graph area background and border.  Note that even if the border is transparent, the
         * width is still taken into account when determining the graph area.
         */
        FillStyle graph_style{White, LineStyle(Black, 1)};
        ///@{
        /// The graph padding inside the border, in surface units.
        double graph_padding_left = 0,
               graph_padding_right = 0,
               graph_padding_top = 2,
               graph_padding_bottom = 2;
        ///@}

        /// The left edge of the legend boxes relative to the right edge of the graph.
        double legend_left = 8;
        /// The top edge of the legend boxes relative to the top edge of the graph.
        double legend_top = 0;
        /// The distance between the right edge of the legend text area and the right edge of the graph.
        double legend_right = 4;
        /// The width of the legend image box, in surface units.
        double legend_box_width = 15;
        /// The height of the legend box, in surface units.
        double legend_box_height = 15;
        /// The gap between the legend box and legend text region, in surface units.
        double legend_box_text_gap = 5;
        /** The maximum height of the box that legend text may occupy (beyond which text is
         * ellipsized), in surface units.  The text will be aligned with the middle of the legend
         * box, and the total legend item height will be the larger of the legend box or text
         * region.  Note that a legend will always contain at least one line of text, and so it is
         * possible for a legend item to exceed this size if this is set smaller than the required
         * height for a line of text.
         */
        double legend_text_max_height = 40;
        /** The font to use for legend text. */
        Pango::FontDescription legend_font = Pango::FontDescription("serif 8");
        /// The vertical space to put between legend items.
        double legend_spacing = 3;

        /// The graph title text, which may contain pango markup.
        std::string title_markup;
        /// The font to use for the title
        Pango::FontDescription title_font = Pango::FontDescription("serif 13");
        /// The space between the top of the image and the top of the title, in surface units.
        double title_padding_top = 2;
        /// The space between the bottom of the title and the top of the graph, in surface units.
        double title_padding_bottom = 2;
        /// The space between the left of the surface and the title, in surface units.
        double title_padding_left = 20;
        /// The space between the right of the surface and the title, in surface units.
        double title_padding_right = 20;

        /** Returns a Cairo::Matrix that transforms coordinates to be graphed into the interior
         * graph area of the image, based on the tmin/tmax/ymin/ymax values given when constructing
         * the Series object.  The returned object will be correct as long as the object's graph_*
         * variables remain unchanged.
         *
         * In the translate space, (tmin,ymin) corresponds to the lower-left corner of the graph
         * space, and (tmax,ymax) corresponds to the top-right corner of the graph space.  (Both
         * "corners" are inside both the border and padding).
         */
        Cairo::Matrix translateGraph() const;

    private:
        Target &target_;
        int tmin_, tmax_;
        double ymin_, ymax_;


        // The buffer to add on the top/bottom edge of the graph (so that lines/points don't run
        // onto the graph outline).  In absolute surface units.
        double graph_buffer_tb_ = 2;
        // The buffer to add to the left/right edges.  Defaults to 0.
        double graph_buffer_lr_ = 0;

        // Whether initializePage() has been called for the current page yet
        bool page_initialized_ = false;

        // Legend items to be re-added when newPage is called
        std::list<std::pair<std::string, LegendPainterCallback_t>> legend_;

        // Vertical offset of the next legend item, from the top of the legend area.
        double legend_next_ = 0;

    protected:
        /** Draws the initial layout (graph box, etc.).  Does nothing if called repeatedly (until
         * the next newPage() call).
         */
        void initializePage();

        /** Clips the given Cairo::Context to the graph space, including the padding and,
         * optionally, the border.  This should be called without a transformation having been
         * applied to the provided context.
         *
         * \param ctx the Cairo::Context for the surface *without* a transformation matrix applied.
         * \param border whether or not the border should be part of the clipped region; if false
         * (the default), the border cannot be drawn on top of; if true, the border can be
         * overwritten.
         */
        void clipToGraph(Cairo::RefPtr<Cairo::Context> ctx, bool border = false) const;

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
