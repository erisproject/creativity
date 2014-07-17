#pragma once
#include <eris/noncopyable.hpp>
#include <eris/Position.hpp>
#include <eris/WrappedPositional.hpp>
#include <gtkmm/drawingarea.h>
#include <mutex>
#include <memory>

namespace creativity {

// forward declarations
class Reader;
class Book;

namespace gui {

class GUI;

/** Gtk drawing area upon which the simulation state visualization is drawn.  */
class GraphArea : public Gtk::DrawingArea, eris::noncopyable {
    public:
        GraphArea() = delete;

        /** Creates a graph area that draws in the rectangle bounded by [`bottom', `top'] on the
         * vertical plane and [`left', `right'] on the horizontal plane.  Specifying a value of
         * bottom or left larger than top or bottom will flip the respective axis.
         */
        GraphArea(const double &top, const double &right, const double &bottom, const double &left, GUI &gui);

        /** Returns a Cairo::Matrix that translates graph coordinates into screen coordinates.  To
         * go the other way, invert this matrix.
         */
        Cairo::Matrix graph_to_canvas() const;

        /** Draws the current set of points/circles, and updates the data in any reader/book
         * dialogs. */
        virtual bool on_draw(const Cairo::RefPtr<Cairo::Context> &cr) override;

        /** The types of points drawn by this class. */
        enum class PointType {
            X, //!< A diagonal cross
            CROSS, //!< A horizontal/vertical cross
            SQUARE //!< A small square
        };

        /** Draws a point.  Drawn points have a fixed size that does not depend on the drawing area
         * size.
         * \param cr the Cairo::Context used for drawing
         * \param trans the transformation matrix that translates from graph space to display space
         * \param x the point x coordinate, in graph space
         * \param y the point y coordinate, in graph space
         * \param type the type of point to draw
         * \param red the red component (0-1) of the colour
         * \param green the green component (0-1) of the colour
         * \param blue the blue component (0-1) of the colour
         * \param alpha the alpha component (0-1) of the colour: 0 is transparent, 1 is opaque.
         * \param scale the scale of the point.  1 (the default) means default size.
         * \param virt whether this is a virtual (i.e. wrapping) point. Internal use only.
         */
        void drawPoint(const Cairo::RefPtr<Cairo::Context> &cr, const Cairo::Matrix &trans, double x, double y,
                const PointType &type, double red, double green, double blue, double alpha, double scale = 1.0, bool virt = false);

        /** Draws a circle.  Circles, however, are in graph space, not drawing area space, so
         * this is actually going to end up drawing ovals (unless the drawing area happens
         * to be perfectly square).
         * \param cr the Cairo::Context used for drawing
         * \param trans the transformation matrix that translates from graph space to display space
         * \param cx the center x coordinate, in graph space
         * \param cy the center y coordinate, in graph space
         * \param r the radius of the circle, in graph space
         * \param red the red component (0-1) of the colour
         * \param green the green component (0-1) of the colour
         * \param blue the blue component (0-1) of the colour
         * \param alpha the alpha component (0-1) of the colour: 0 is transparent, 1 is opaque.
         */
        void drawCircle(const Cairo::RefPtr<Cairo::Context> &cr, const Cairo::Matrix &trans,
                double cx, double cy, double r,
                double red, double green, double blue, double alpha);

        /// Spacing of axes tick marks.  2.0 means ticks at 2, 4, 6, etc.
        double tick_space = 1.0;
        /// The size of tick marks, in pixels.
        double tick_size = 6.0;
        /// Every `tick_big`th tick will be triple-sized
        int tick_big = 5;

        /** The radius of markers representing points, in pixels.  For example, a value of 5 means
         * that CROSS points will have horizontal and vertical lines extending a distance of 5,
         * while X points will have diagonal lines with a length of 5.  For SQUARE points the
         * diagonals of the square will have length 10 (and so edges will be \f$\5 \sqrt{2}\f$
         * pixels long).
         */
        double point_size = 5.0;
    private:
        // The parent GUI
        GUI &gui_;
        // Helper object for doing wrapping calculations
        eris::WrappedPositionalBase wpb_;
        // The bounds of the graph
        struct {
            double top; ///< the top graph boundary 
            double right; ///< the right graph boundary
            double bottom; ///< the bottom graph boundary
            double left; ///< the left graph boundary
        } bounds_;

        // Cache of image surfaces; these save redrawing when transitioning between states
        std::vector<Cairo::RefPtr<Cairo::ImageSurface>> drawing_cache_;
        int drawing_cache_width_ = 0, drawing_cache_height_ = 0;

        // Sets up one or more (wrapping) lines between two points, but does not actually draw it
        // with stroke().
        void drawWrappingLine(const Cairo::RefPtr<Cairo::Context> &cr, const Cairo::Matrix &trans, const eris::Position &from, const eris::Position &to);
};

} }
