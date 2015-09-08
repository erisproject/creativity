#pragma once
#include <eris/noncopyable.hpp>
#include <eris/WrappedPositional.hpp>
#include <gtkmm/drawingarea.h>
#include <cairomm/matrix.h>
#include <cairomm/pattern.h>
#include <cairomm/refptr.h>
#include <cairomm/surface.h>
#include <cstddef>
#include <vector>

namespace Cairo { class Context; }

namespace creativity { namespace gui {

class GUI;

/** Gtk drawing area upon which the simulation state visualization is drawn.  */
class GraphArea : public Gtk::DrawingArea, eris::noncopyable {
    public:
        GraphArea() = delete;

        /** Creates a graph area that draws according to the boundary of the given GUI object.
         */
        GraphArea(GUI &gui);

        /** Returns a Cairo::Matrix that translates graph coordinates into canvas (i.e. on-screen
         * pixel) coordinates.  To go the other way, call canvas_to_graph(), or .invert() the
         * returned matrix.
         *
         * Note that the graph is always 2-dimensional; when the simulation is 2-dimensional,
         * graph positions are exactly simulation positions.  When 1-dimensional, positions are
         * mapped from 1D into a graph circle.  When 3-dimensional, the graph position is the
         * simulation position in the `design.subdimensionX` and `design.subdimensionY` indices.
         */
        Cairo::Matrix graph_to_canvas() const;

        /** Returns a Cairo::Matrix that tranlates canvas coordinates into graph coordinates.  To go
         * the other way, call graph_to_canvas(), or .invert() the returned matrix.
         */
        Cairo::Matrix canvas_to_graph() const;

        /** Returns the graph position (with each dimension between -boundary and +boundary) from a
         * simulation Position.
         *
         * If the simulation is two dimensional, the graph position equals the simulation Position.
         * If one-dimensional, the graph position is the position on a circle of radius 0.9*boundary
         * 1 dimension, where position 0 is the right-most point of the circle, boundary/2 is the
         * top, boundary (and -boundary) is the left-most point, etc.
         *
         * If 3- (or more) dimensional, the `design.subdimensionX` and `design.subdimensionY`
         * dimension indices are used.
         */
        eris::Position graph_position(const eris::Position &p);

        /** Draws the current set of points/circles, and updates the data in any reader/book
         * dialogs. */
        virtual bool on_draw(const Cairo::RefPtr<Cairo::Context> &cr) override;

        /** The types of points drawn by this class. */
        enum class PointType {
            X, //!< A diagonal cross, used for readers.
            CROSS, //!< A horizontal/vertical cross, used for books.
            SQUARE //!< A small square.  Currently unused.
        };

        /// Typedef of Colour to the appropriate ref-pointed Cairo class
        using Colour = Cairo::RefPtr<Cairo::SolidPattern>;

        /** Draws a point.  Drawn points have a fixed size that does not depend on the drawing area
         * size.
         * \param cr the Cairo::Context used for drawing
         * \param trans the transformation matrix that translates from graph space to display space
         * \param point the location of the point, in simulation space
         * \param type the type of point to draw
         * \param colour the colour (and alpha channel) of the point
         * \param radius the radius of the point.  1 (the default) means default size.
         * \param line_width the line width to use to draw the point indicator
         */
        void drawPoint(const Cairo::RefPtr<Cairo::Context> &cr, const Cairo::Matrix &trans, const eris::Position &point,
                const PointType &type, const Colour &colour, double radius, double line_width);

        /** Draws a circle in graph space, not canvas space, so this is actually going to end up
         * drawing ovals (unless the canvas happens to be perfectly square).  The circle can also
         * have a radial line drawn at a random angle, if given a non-zero radial stroke width.  The
         * angle is calculated from the given coordinates and radius, and so will not change across
         * redraws.
         *
         * \param cr the Cairo::Context used for drawing
         * \param trans the transformation matrix that translates from graph space to canvas space
         * \param cx the centre x coordinate, in graph space
         * \param cy the centre y coordinate, in graph space
         * \param r the radius of the circle, in graph space
         * \param colour the colour (and possibly alpha value) of the circle.
         * \param stroke_width the width of the stroke used to draw the circle.
         * \param radial_stroke_width the width of the stroke used to draw a line from the centre to
         * the circle circumference (at a random angle).  If 0, no radial is drawn.
         *
         * Also note that circles currently do *not* wrap, unlike points and lines.
         *
         * \sa drawCanvasCircle for drawing a circle in canvas space, i.e. something that will
         * always show up to the user as a circle.
         */
        void drawGraphCircle(const Cairo::RefPtr<Cairo::Context> &cr, const Cairo::Matrix &trans,
                double cx, double cy, double r, const Colour &colour, double stroke_width, double radial_stroke_width = 0.0);

        /** Draws a circle in canvas space, which will always appear circular on the GUI element
         * (but won't be circular in graph space).
         *
         * A radial line will also be drawn (if `radial_stroke_width > 0`) from the center of the
         * circle to a random point of the circumference.  The pseudo-random angle is calculated
         * using the given center point and radius so that redrawing and resizing the canvas will
         * not change the angle.
         *
         * \param cr the Cairo::Context used for drawing
         * \param trans the transformation matrix that transforms graph points to canvas points
         * \param centre the centre location, in simulation (not graph or canvas!) space
         * \param r the radius of the circle, in canvas (not graph!) space
         * \param colour the colour (and possibly alpha value) of the circle.
         * \param stroke_width the width of the stroke used to draw the circle
         * \param radial_stroke_width the width of the stroke used to draw a line from the centre to
         * the circle circumference (at a random angle).  If 0, no radial is drawn.
         *
         * \sa drawCanvasCircle for drawing a circle in canvas space (which generally looks like an
         * oval on the screen, unless the GraphArea canvas happens to be exactly square).
         */
        void drawCanvasCircle(const Cairo::RefPtr<Cairo::Context> &cr, const Cairo::Matrix &trans,
                const eris::Position &centre, double r, const Colour &colour, double stroke_width, double radial_stroke_width = 0.0);

        /** Draws the shortest line between two points, allowing the line to wrap across simulation
         * boundaries.  Note that when a line wraps this actually involves drawing several line
         * segments.
         *
         * \param cr the Cairo::Context on which to draw
         * \param trans the transformation matrix that converts graph points to canvas points
         * \param from the start position of the line, in simulation space
         * \param to the end position of the line, in simulation space
         * \param min_length the minimum length line to draw; if the on-screen line length (in
         * pixels) is less than this, it will not be drawn.  If omitted, defaults to 0 (no minimum).
         * \param start_colour the start colour of the line (if not specified, the stroke colour is
         * not changed), for creating a line with a gradient.
         * \param end_colour the ending colour of the line, for creating a line with a gradient.
         */
        void drawLine(const Cairo::RefPtr<Cairo::Context> &cr, const Cairo::Matrix &trans,
                const eris::Position &from, const eris::Position &to,
                double min_length = 0.0, Colour start_colour = Colour(), Colour end_colour = Colour());

        /** Draws a pair of graph axes going to `-boundary` to `+boundary` in each dimension.
         *
         * \param cr the Cairo::Context on which to draw
         * \param trans the transformation matrix that converts graph points to canvas points
         */
        void drawAxes(const Cairo::RefPtr<Cairo::Context> &cr, const Cairo::Matrix &trans);

        /** Draws a circular single axis for a 1D model that goes from 0 (at the right) rotating
         * counterclockwise through to `boundary` and clockwise to `-boundary` (both are the
         * leftmost point).
         *
         * \param cr the Cairo::Context on which to draw
         * \param trans the transformation matrix that converts graph points to canvas points
         */
        void drawCircularAxis(const Cairo::RefPtr<Cairo::Context> &cr, const Cairo::Matrix &trans);

        /** Struct containing various graph properties such as colours, line widths, and dash
         * styles.
         */
        struct {
            /// Struct controlling whether various elements are drawn in the graph
            struct {
                bool
                    reader = true, ///< Whether to draw agents
                    friendship = true, ///< Whether to draw friendship links
                    movement = true, ///< Whether to draw a reader movement trail
                    book_live = true, ///< Whether to draw on-market books
                    book_dead = true, ///< Whether to draw off-market books
                    book_public = true, ///< Whether to draw public-market books
                    author_live = true, ///< Whether to draw on-market authorship lines
                    author_dead = true, ///< Whether to draw off-market authorship lines
                    author_public = true, ///< Whether to draw public-market authorship lines
                    reading = true, ///< Whether to draw lines for newly obtained books
                    utility_gain = true, ///< Whether to draw utility gain circles
                    utility_loss = true, ///< Whether to draw utility loss circles
                    axes = true; ///< Whether to draw the graph axes
            } enabled;

            /// The colours of various visualization elements
            struct {
                Colour
                    /// The colour of private, on-market books
                    book_live = Cairo::SolidPattern::create_rgb(0, .4, 1),
                    /// The colour of off-market books
                    book_dead = Cairo::SolidPattern::create_rgb(0.5, 0.5, 0.5),
                    /// The colour of public-market books
                    book_public = Cairo::SolidPattern::create_rgb(0, 1, .4),
                    /// Authorship lines from author to book for on-market books
                    author_live = Cairo::SolidPattern::create_rgba(0.5, 0.2, 0.5, 0.5),
                    /// Authorship lines from author to book for off-market books
                    author_dead = Cairo::SolidPattern::create_rgba(0.75, 0.5, 0.75, 0.5),
                    /// Authorship lines from author to book for public-market books
                    author_public = Cairo::SolidPattern::create_rgba(0.5, 0.2, 0.5, 0.25),
                    /// The colour of lines from readers to newly purchased books
                    reading = Cairo::SolidPattern::create_rgba(1, 0.55, 0, 0.5),
                    /// The colour of agents (readers/authors)
                    reader = Cairo::SolidPattern::create_rgb(1, 0, 0),
                    /// The colour of friendship links between readers
                    friendship = Cairo::SolidPattern::create_rgba(0.75, 0, 0, 0.5),
                    /** The colour of the movement trail; the actual trail will be a gradient from
                     * fully transparent to this colour.
                     */
                    movement = Cairo::SolidPattern::create_rgba(1, 0, 0, 0.75),
                    /// The colour of the utility circle (for readers with utility > 1000)
                    utility_gain = Cairo::SolidPattern::create_rgba(.133, .545, .133, 0.5),
                    /// The colour of the utility loss circle (for readers with utility < 1000)
                    utility_loss = Cairo::SolidPattern::create_rgba(.545, .133, .133, 0.5),
                    /// The colour of the graph axes and tick marks
                    axes = Cairo::SolidPattern::create_rgb(0, 0, 0),
                    /// The colour of the background
                    background = Cairo::SolidPattern::create_rgb(1, 1, 1);
            } colour;

            /** The dash pattern of various visualization elements.  An empty vector results in a
             * solid line; otherwise vector elements are alternating lengths of on and off line
             * segments along the path.
             */
            struct {
                std::vector<double>
                    /// The dash pattern of author lines to on-market books
                    author_live,
                    /// The dash pattern of author lines to off-market books
                    author_dead{{3.0, 2.0}},
                    /// The dash pattern of author lines to public-market books
                    author_public{{3.0, 2.0}},
                    /// The dash pattern of lines from readers to newly purchased books
                    reading,
                    /// The dash pattern of friendship links between readers
                    friendship{{12.0, 3.0}},
                    /// The dash pattern of reader movement trails
                    movement,
                    /// The dash pattern of the utility circle (for readers with utility > 1000)
                    utility;
            } dash;

            /// The stroke width used to draw various visualization elements
            struct {
                double
                    book_live = 2.0, ///< Stroke width of live books
                    book_dead = 2.0, ///< Stroke width of dead books
                    book_public = 2.0, ///< Stroke width of public books
                    author_live = 2.0, ///< Stroke width of author-to-live book lines
                    author_dead = 2.0, ///< Stroke width of author-to-dead book lines
                    author_public = 2.0, ///< Stroke width of author-to-dead book lines
                    reading = 2.0, ///< Width of lines from readers to newly-obtained books
                    reader = 2.0, ///< Stroke width of agent icons
                    friendship = 1.0, ///< Stroke width of reader-to-reader friendship links
                    movement = 3.0, ///< Stroke width of movement path
                    utility = 2.0, ///< Stroke width of utility circles
                    utility_radial = 1.0, ///< Stroke width of the radial line of utility circles
                    axes = 2.0, ///< Stroke width of the graph axes
                    axes_ticks = 1.0, ///< Stroke width of graph axes tick marks
                    axes_ticks_big = 2.0; ///< Stroke width of every nth tick mark
            } stroke_width;

            /// The sizes of visual elements
            struct {
                double
                    reader = 7.071, ///< Radius of a reader "x" symbol
                    /** Minimum radius of a book "+" symbol.  Note that newer books are scaled to a
                     * multiple of this size.
                     */
                    book = 5.0,
                    /** Scaling parameter 'a' for live books; books are scaled to the above
                     * book_live parameter times \f$\max(1.0, a - b \times age)\f$.
                     */
                    book_scale_a = 3.0,
                    /** Scaling parameter 'b' for live books; books are scaled to the above
                     * book_live parameter times \f$\max(1.0, a - b \times age)\f$.
                     */
                    book_scale_b = 0.3,
                    /** The radius multiple for the reader utility gain circle; the circle is drawn
                     * with a radius around the reader of \f$r*ln(u-999)\f$, where this value is
                     * \f$r\f$, but only for readers with utility above 1000.
                     */
                    utility_gain_scale = 5.0,
                    /** The radius multiple for the reader utility loss circle; the circle is drawn
                     * with a radius around the reader of \f$r*ln(1001-u)\f$, where this value is
                     * \f$r\f$, but only for readers with utility below 1000.
                     */
                    utility_loss_scale = 5.0,
                    /** Length of regular tick marks, in pixels.  The tick extends half this
                     * distance away from the axis in both directions.
                     */
                    tick = 6.0,
                    /// Length of big tick marks (every 5th tick, by default), in pixels.
                    tick_big = 18.0;

            } size;

            /// The visual design of types of points (readers and books)
            struct {
                PointType
                    /// The point type of readers
                    reader = PointType::X,
                    /// The point type of live books
                    book_live = PointType::CROSS,
                    /// The point type of dead books
                    book_dead = PointType::CROSS,
                    /// The point type of public books
                    book_public = PointType::CROSS;

                double
                    /// Spacing of axes tick marks.  2.0 means ticks at 2, 4, 6, etc.
                    tick_every = 1.0,
                    /// Alpha multiplier for the starting colour of the movement line gradient.
                    movement_alpha_multiplier = 0.1;

                unsigned int
                    /// Every `style.tick_big`th tick will be `size.tick_big`-sized
                    tick_big = 5;
            } style;

            /** If displaying a 3 (or more) dimensional simulation, this is the dimension index to
             * display on the X axis.
             */
            size_t subdimensionX = 0;
            /** If displaying a 3 (or more) dimensional simulation, this is the dimension index to
             * display on the Y axis.
             */
            size_t subdimensionY = 1;

        } design;

        /** Resets the drawn image cache.  This should be called after changing anything in
         * `.design`.  This does not, however, trigger a redraw, so this should typically be
         * followed by a call to `.queue_draw()`.
         */
        void resetCache();

    protected:
        /** Called by drawPoint() with graph coordinates to draw individual canvas points.  This
         * will be called from 1 to 4 times per drawPoint() call: once for the actual point, and
         * potential other calls for the wrapped version when the point is close to an edge or
         * corner.
         *
         * The given `x` and `y` values are the graph coordinates (which will be translated in
         * canvas coordinates).
         */
        void drawPointSingle(const Cairo::RefPtr<Cairo::Context> &cr, const Cairo::Matrix &trans,
                double x, double y, const PointType &type, const Colour &colour, double radius, double line_width);

        /// The parent GUI
        GUI &gui_;
        /// Helper object for doing wrapping calculations
        eris::WrappedPositionalBase wpb_;

        /// Cache of image surfaces; these save redrawing when transitioning between states
        std::vector<Cairo::RefPtr<Cairo::ImageSurface>> drawing_cache_;
        /// The width of the elements in drawing_cache_.
        int drawing_cache_width_ = -1;
        /// The height of the elements in drawing_cache_.
        int drawing_cache_height_ = -1;

};

} }
