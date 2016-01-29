#pragma once
#include "creativity/data/graph/style.hpp"
#include "creativity/data/graph/Target.hpp"
#include <cairomm/context.h>
#include <pangomm/fontdescription.h>
#include <pangomm/font.h>
#include <map>
#include <set>

#ifndef DOXYGEN_SHOULD_SEE_THIS
#if CAIROMM_MAJOR_VERSION < 1 || (CAIROMM_MAJOR_VERSION == 1 && CAIROMM_MINOR_VERSION <= 12)
namespace Cairo {

/** RAII-style context save/restore class.  context->save() is called
 * automatically when the object is created, and context->restore() is called
 * when the object is destroyed.  This allows you to write code such as:
 *
 *     // context initial state
 *     {
 *         Cairo::SaveGuard saver(context);
 *         ... // manipulate context
 *     }
 *     // context is restored to initial state
 */
class SaveGuard final {
  public:
    /// Constructor: the context is saved
    explicit SaveGuard(RefPtr<Context> context) : ctx_{context} { ctx_->save(); }
    /// Copy constructor deleted
    SaveGuard(const SaveGuard &) = delete;
    /// Destructor; the context is restored
    ~SaveGuard() { ctx_->restore(); }
  private:
    RefPtr<Context> ctx_;
};
}
#endif
#endif

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

        /// Destructor; calls finishPage()
        virtual ~Series();

        /** Writes everything added to the graph to the current page.  Once this method is called,
         * no more graph elements may be added until newPage() is called to begin a new graph.
         *
         * This method is called implicitly by `newPage()` and by the destructor; as such it is
         * rarely required to call it directly.
         */
        void finishPage();

        /** Writes the current graph to disk by calling finishPage(), if necessary, and begins a new
         * graph.  If a multi-page surface format, this will be a new page in the output file;
         * otherwise it will be a new file.
         */
        void newPage();

        /// The default line style for addLine(): a 1-unit wide, solid black line.
        static const LineStyle default_line_style;

        /** Plots the pairs of coordinates in the given map to the plot.  Keys should be the `t`
         * values.  If any sequential t values are missing or have non-finite values, the line will
         * be broken around the missing values.
         *
         * \param points the map of time to value coordinates to plot
         * \param style the LineStyle; defaults to a 1-unit-wide black line.
         */
        void addLine(const std::map<int, double> &points,
                const LineStyle &style = default_line_style);

        /** The default fill style for a region added with addRegion(): a 1/4 opacity solid blue
         * background, with no border.
         */
        static const FillStyle default_region_style;

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
        void addRegion(const std::map<int, std::pair<double, double>> &intervals,
                const FillStyle &style = default_region_style
                );

        /** The default legend box background and/or border style: white background with a 1-unit
         * wide, solid black border.
         */
        static const RectangleStyle default_legend_box_style;

        /** Adds a simple legend item.  The legend box has two main drawing components; each can be
         * transparent to suppress the element.  If the page has not yet been initialized *and* the
         * `preserve` value is true, the legend item will not be drawn immediately but rather will
         * wait until the current page is initialized.  Otherwise, the legend item will be drawn
         * immediately, initializing the page first if required.
         *
         * \param markup the text (with pango markup) to display for the legend item
         * \param lower fill style for the bottom half of the box.  The border component of the fill
         * style is drawn across the middle of the box.  To specify just a line across the middle,
         * simple provide a LineStyle for this argument (the implicit conversion results in a
         * transparent fill and the given line style).
         * \param preserve whether this legend item should be automatically added to new pages when
         * newPage() is called.
         * \param style the fill and border style for the legend box.  Defaults to a transparent
         * background with a 1-unit-wide black border.  Note that the background fill is always
         * contained entirely within the border, even if the border is transparent.
         *
         * The background fill and border is drawn first, followed by the lower area.
         */
        void addLegendItem(std::string markup,
                FillStyle lower,
                bool preserve = true,
                RectangleStyle style = default_legend_box_style);

        /// Alias for the function callback for custom legend support
        using LegendPainterCallback_t = std::function<void(Cairo::RefPtr<Cairo::Context>, const Cairo::Matrix&)>;

        /** Adds a legend item with custom box drawing.  The custom callback function is called with
         * a Cairo::Context object to allow the caller to draw the box as desired.  If the current
         * page is uninitialized *and* `preserve` is true, the legend item drawing will be delayed
         * until page initialization; otherwise the legend item is drawn immediately, initializing
         * the page if required.
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
         * \param style a RectangleStyle from which to read the border style and background fill colour
         * \param clip if true, clip the context to the interior of the rectangle before returning
         */
        static void drawRectangle(Cairo::RefPtr<Cairo::Context> ctx, double width, double height, const RectangleStyle &style, bool clip = false);

        /// The background colour of the page
        RGBA background_colour = White;

        /** The line style (including length) of tick marks.  Defaults to a black, 0.5-width,
         * 3-length line.  The `length` parameter of the style controls the length of tick marks; if
         * 0, tick marks are suppressed. */
        LineStyle tick_style = LineStyle(Black, 0.5, 3);

        /** The space between the tick line and the corresponding tick value text.  Only applies
         * when both are visible (i.e. non-zero `tick_style` width and length, and non-transparent
         * `tick_font_colour`).
         */
        double tick_label_space = 2;

        /** The line style to draw for ticks grid lines on the graph itself (*not* the tick marks
         * themselves, which show up outside the graph area). Defaults to transparent (i.e. no grid
         * lines). */
        LineStyle tick_grid_style = LineStyle(Transparent, 0.5);

        /** The font colour for tick values.  If completely transparent, tick values are not drawn. */
        RGBA tick_font_colour = Black;

        /** The font for tick values */
        Pango::FontDescription tick_font = Pango::FontDescription("serif 6");

        /** The graph area background and border.  Note that even if the border is transparent, the
         * width is still taken into account when determining the graph area.
         */
        RectangleStyle graph_style{White, LineStyle(Black, 1)};
        ///@{
        /// The graph padding inside the border, in surface units.
        double graph_padding_left = 0,
               graph_padding_right = 0,
               graph_padding_top = 2,
               graph_padding_bottom = 2;
        ///@}

        ///@{
        /// The minimum margin between components of the page and the edge of the page on each side.
        double margin_top = 5,
               margin_right = 5,
               margin_bottom = 5,
               margin_left = 5;
        ///@}

        /** The possible positions of the legend area (before adjustment by legend_graph_space,
         * legend_x_offset and legend_y_offset).  Note that the pre-offset position is such that the
         * legend border lies exactly adjacent to graph border (you can set `legend_graph_space` to
         * `-graph_style.border.thickness` to make the borders overlap).
         */
        enum class LegendPosition {
            None, ///< The legend is suppressed entirely
            OutsideTop, ///< To the right of the graph, aligned relative to the graph top
            OutsideMiddle, ///< To the right of the graph, centered relative to the graph
            OutsideBottom, ///< To the right of the graph, aligned relative to the graph bottom
            TopRight, ///< Inside the graph, in the top-right corner
            MiddleRight, ///< Inside the graph, in the middle of the right side
            BottomRight, ///< Inside the graph, in the bottom-right corner
            TopCenter, ///< Inside the graph, in the middle of the top edge
            MiddleCenter, ///< Inside the graph, in the middle of the graph
            BottomCenter, ///< Inside the graph, in the middle of the bottom edge
            TopLeft, ///< Inside the graph, in the top-left corner
            MiddleLeft, ///< Inside the graph, in the middle of the left-hand side
            BottomLeft ///< Inside the graph, in the bottom-left corner
        };

        /// The position in which to draw the legend.
        LegendPosition legend_position = LegendPosition::OutsideTop;

        /** The size of the gap between the legend area and the graph side(s) with which the legend
         * is aligned.  The precise adjustment depends on the LegendPosition:
         *
         * - OutsideTop, OutsideMiddle, OutsideBottom -- larger values move the legend right
         * - TopRight -- larger values move the legend left and down
         * - MiddleRight -- larger values move the legend left
         * - BottomRight -- larger values move the legend left and up
         * - TopCenter -- larger values move the legend down
         * - MiddleCenter -- no effect
         * - BottomCenter -- larger values move the legend up
         * - TopLeft -- larger values move the legend down and right
         * - MiddleLeft -- larger values move the legend right
         * - BottomLeft -- larger values move the legend right and up
         */
        double legend_graph_space = 5.0;
        ///@{
        /** After adjusting the position according to the `legend_graph_space` value, the x and y
         * coordinates are also adjusted by these values.  Larger values are further right (x) or
         * down (y).  Note that any offset specified here is not taken into account when calculating
         * the page layout: large positive or negative values can move the legend into the page
         * margin or off the page entirely.
         */
        double legend_x_offset = 0;
        double legend_y_offset = 0;
        ///@}
        /// The width of the legend image box, in surface units.
        double legend_box_width = 15;
        /// The height of the legend box, in surface units.
        double legend_box_height = 15;
        /// The gap between the legend box and legend text region, in surface units.
        double legend_box_text_gap = 5;
        ///@{
        /// The legend box area padding inside the border, in surface units.
        double legend_padding_left = 3,
               legend_padding_right = 3,
               legend_padding_top = 3,
               legend_padding_bottom = 3;
        ///@}
        /// The style for a box encompassing the entire legend area.  Default is a white background
        /// and 1-unit wide black border.
        RectangleStyle legend_style{White, Black};
        /** The maximum width of legend text, in surface units.  The actual width can be less than
         * this, depending on the size of the text and taking into account any wrapping.  If text is
         * too wide (for example, one unwrappable line exceeds this value) it will be cut off; if
         * the text is too long overall, it will be ellipsized.
         */
        double legend_text_max_width = 60;
        /** The maximum height of the box that legend text may occupy (beyond which text is
         * ellipsized), in surface units.  The text will be aligned with the middle of the legend
         * box, and the total legend item height will be the larger of the legend box or text
         * region.  Note that a legend will always contain at least one line of text, and so it is
         * possible for a legend item to exceed this size if this is set smaller than the required
         * height for a line of text.
         */
        double legend_text_max_height = 40;
        /// The vertical space to put between legend items.
        double legend_spacing = 3;
        /** The font to use for legend text. */
        Pango::FontDescription legend_font = Pango::FontDescription("serif 6");

        /** The graph title text, which may contain pango markup.  Will be drawn in opaque black,
         * but markup in the string can be used to change this.
         */
        std::string title_markup;
        /// The font to use for the title
        Pango::FontDescription title_font = Pango::FontDescription("serif 12");
        /** The space between the bottom of the title and the top of the graph, in surface units.
         * Only applies when title_markup is not an empty string.  To adjust the space between the
         * title and the top of the page, adjust `margin_top`.
         */
        double title_padding_bottom = 2;
        /// The minimum space between the left of the surface and the title, in surface units.
        double title_padding_left = 20;
        /// The minimum space between the right of the surface and the title, in surface units.
        double title_padding_right = 20;

        /** The thickness of the region drawn for "stray" region observations, that is, observations
         * that are neither preceeded nor followed by another finite region observation.
         */
        double stray_thickness = 1;

        /** X-axis (t) tick mark locations, calculated during construction.  Can be recalculated by
         * calling recalcTicks().  If cleared, no horizontal axis ticks are displayed.
         */
        std::set<int> t_ticks;

        /** Y-axis tick mark locations, calculated during construction.  Can be recalculated by
         * calling recalcTicks().  If cleared, no vertical axis ticks are displayed.
         */
        std::set<double> y_ticks;

        /** X-axis (t) grid mark locations, calculated during construction.  Can be recalculated by
         * calling recalcTicks().  If cleared, no vertical grid lines are drawn.
         *
         * This will, by default, be equal to t_ticks if recalcTicks() is called with TickEnds::Add;
         * for the other end modes, the grid end points can differ.
         */
        std::set<int> t_grid;

        /** Y-axis grid mark locations, calculated during construction.  Can be recalculated by
         * calling recalcTicks().  If cleared, no horizontal grid lines are drawn.
         *
         * This will, by default, be equal to y_ticks if recalcTicks() is called with TickEnds::None;
         * for the other end modes, the grid end points can differ.
         */
        std::set<int> y_grid;

        /** Different modes for what to do with the tmin/tmax/ymin/ymax graph end points given
         * during construction. */
        enum class TickEnds {
            /** If an end point isn't equal to or very close to a calculated tick mark, add it.
             * This can result in up to two more than the requested maximum number of ticks.  Very
             * close is considered to be within one-quarter of a tick from the first or last
             * calculated tick.
             *
             * For example, if the end points result in calculated tick marks are 5 10 15 20 25 30
             * 35, with TickEnds::Add an extra tick would be added to the beginning for a minimum
             * value in (0, 3.75], and would be added to the end for a maximum value in [36.25, 40);
             * ticks are not added for end points within 1.25 (1/4 of a tick) of the calculated
             * minimum/maximum ticks.
             */
            Add,
            Replace, ///< The first and last calculated tick marks are removed and replaced with the min/max graph axis values.
            None ///< Don't do anything with tick end values.  Value range ends will only show up if the ends happen to coincide with multiple of the tick increments
        };

        /** Calculates where tick marks and grid lines should go.  This will be applied the next
         * time the graph is drawn, and so has to be called before adding any graph elements to have
         * an effect on the current page.
         *
         * This method simply calls recalcTicksT() and recalcTicksY().
         *
         * \param xmax the maximum number of tick marks on the x axis.  Defaults to 10
         * \param ymax the maximum number of tick marks on the y axis.  Defaults to 8
         * \param end_mode what to do with the ends (i.e. tmin, tmax, ymin, ymax given during
         * construction).  The default is TickEnds::None.  Grid lines are not affected by this
         * option.
         */
        void recalcTicks(unsigned xmax = 10, unsigned ymax = 8, TickEnds end_mode = TickEnds::None);

        /** Recalculates just the t tick marks and grid lines on the graph, leaving y ticks
         * unchanged.  The tick marks are calculated to have a round increment which is the smallest
         * of 1, 2, 4, 5, 10, 20, 25, 40, 50, 100, 200, 250, 400, 500, ... such that no more than
         * `max` ticks are needed to cover the range of t data.
         *
         * \param max the maximum number of tick marks on the t/x axis.  Defaults to 10.
         * \param end_mode what to do with the ends (i.e. tmin, tmax) for tick marks.  The default
         * is TickEnds::None.  Grid lines are not affected by this option.
         */
        void recalcTicksT(unsigned max = 10, TickEnds end_mode = TickEnds::None);

        /** Recalculates just the y tick marks and grid lines on the graph, leaving t ticks
         * unchanged.  The space between tick marks is calculated to be smallest value of (1, 2,
         * 2.5, 4, 5) times a power of ten that results in no more than `max` ticks to cover the
         * span of the data given during construction.
         *
         * \param max the maximum number of tick marks on the y axis.  Defaults to 8.
         * \param end_mode what to do with the ends (i.e. ymin, ymax).  The default is
         * TickEnds::None.  Grid lines are not affected by this option.
         */
        void recalcTicksY(unsigned max = 8, TickEnds end_mode = TickEnds::None);

        /** Updates the t and y graph extents.  This must be called before the current page is
         * initialized--that is, either after construction, or after a newPage() call before any
         * items have been added to the graph.
         *
         * \param tmin the new minimum t value
         * \param tmax the new maximum t value
         * \param ymin the new minimum y value
         * \param ymax the new maximum y value
         * \param retick if true or omitted, call recalcTicks() after updating the extents (with its
         * default arguments).  Set to false if using custom tick calculations.
         * \throws std::logic_error if called after the current page has been initialized
         * \throws std::invalid_argument if called with invalid extents (tmin not less than tmax, or
         * ymin not less than ymax).
         */
        void setExtents(int tmin, int tmax, double ymin, double ymax, bool retick = true);

    private:
        Target &target_;
        int tmin_, tmax_;
        double ymin_, ymax_;

        // If this is true, the current page is done.
        bool page_finished_ = false;

        using LegendItem_t = std::tuple<bool, std::string, LegendPainterCallback_t>;
        // Legend items to add to the page; <preserve, title, callback> tuples.
        std::list<LegendItem_t> legend_items_;
        // Draw callbacks for the graph items; each will be called with the ctx already clipped to
        // the graph region, and surrounded with a save/restore pair.
        std::list<std::function<void(Cairo::RefPtr<Cairo::Context> ctx, const Cairo::Matrix&)>> draw_;

        // Implementation callback for addLegendItem(markup, lower, preserve, bg).
        void drawSimpleLegendItem(Cairo::RefPtr<Cairo::Context> ctx, const Cairo::Matrix &to_box, FillStyle lower, RectangleStyle bg) const;

        // Draws a legend item at the current position of the context, returning its size (width and height).  If dry_run is true, just returns the size.
        std::pair<double, double> drawLegendItem(Cairo::RefPtr<Cairo::Context> ctx, const LegendItem_t &item, bool dry_run = false) const;

        // Implementation method of addLine() that does the actual drawing
        void drawLine(Cairo::RefPtr<Cairo::Context> ctx, const Cairo::Matrix &translate_graph, const std::map<int, double> &points, const LineStyle &style) const;

        // Implementation method of addRegion() that does the actual drawing
        void drawRegion(Cairo::RefPtr<Cairo::Context> ctx, const Cairo::Matrix &translate_graph, const std::map<int, std::pair<double,double>> &points, const FillStyle &style) const;

    protected:
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
         * \param strays a map to store a "stray" (that is, single-time-period) observation in (only
         * when first == last).
         */
        static void boundRegion(
                Cairo::RefPtr<Cairo::Context> ctx,
                std::map<int, std::pair<double, double>>::const_iterator first,
                std::map<int, std::pair<double, double>>::const_iterator last,
                std::map<int, std::pair<double, double>> &strays);
};

}}}
