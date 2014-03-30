#pragma once
#include <eris/Simulation.hpp>
#include <eris/noncopyable.hpp>
#include <gtkmm/drawingarea.h>
#include <mutex>
#include <memory>

namespace creativity {

// forward declarations
class GUI;

class GUIGraphArea : public Gtk::DrawingArea, eris::noncopyable {
    public:
        /** Creates a graph area that draws in the rectangle bounded by [`bottom', `top'] on the
         * vertical plane and [`left', `right'] on the horizontal plane.  Specifying a value of
         * bottom or left larger than top or bottom will flip the respective axis.
         */
        GUIGraphArea(const double &top, const double &right, const double &bottom, const double &left,
                std::shared_ptr<eris::Simulation> sim, GUI &gui);

        /** Returns a Cairo::Matrix that translates graph coordinates into screen coordinates.  To
         * go the other way, invert this matrix.
         */
        Cairo::Matrix graph_to_canvas() const;

        /** Draws the current set of points/circles. */
        virtual bool on_draw(const Cairo::RefPtr<Cairo::Context> &cr) override;

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
         * \param scale the scale of the point.  1 (the default) means default size.
         */
        void drawPoint(const Cairo::RefPtr<Cairo::Context> &cr, const Cairo::Matrix &trans,
                double x, double y, const PointType &type, double scale = 1.0);

        /** Circle types supported by addCircle() */
        enum class CircleType {
            A, //!< Circle style A
            B //!< Circle style B, will look different from A
        };

        /** Draws a circle.  Circles, however, are in graph space, not drawing area space, so
         * this is actually going to end up drawing ovals (unless the drawing area happens
         * to be perfectly square).
         * \param cr the Cairo::Context used for drawing
         * \param trans the transformation matrix that translates from graph space to display space
         * \param cx the center x coordinate, in graph space
         * \param cy the center y coordinate, in graph space
         * \param r the radius of the circle, in graph space
         * \param type the type of circle to draw
         */
        void drawCircle(const Cairo::RefPtr<Cairo::Context> &cr, const Cairo::Matrix &trans,
                const double &cx, const double &cy, const double &r, const CircleType &type);

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
        // simulation object
        std::shared_ptr<eris::Simulation> sim_;
        // The parent GUI
        GUI &gui_;
        // The bounds of the graph
        std::array<double, 4> bounds_;
        // Constants for access into bounds_
        static constexpr size_t TOP = 0, RIGHT = 1, BOTTOM = 2, LEFT = 3;
};

}
