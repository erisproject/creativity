#pragma once
#include <eris/Simulation.hpp>
#include <eris/noncopyable.hpp>
#include <gtkmm/drawingarea.h>
#include <mutex>
#include <memory>

namespace creativity {

class GUI;

class GUIGraphArea : public Gtk::DrawingArea, eris::noncopyable {
    public:
        /** Creates a graph area that draws in the rectangle bounded by [`bottom', `top'] on the
         * vertical plane and [`left', `right'] on the horizontal plane.  Specifying a value of
         * bottom or left larger than top or bottom will flip the respective axis.
         */
        GUIGraphArea(const double &top, const double &right, const double &bottom, const double &left,
                GUI &gui, std::shared_ptr<eris::Simulation> sim);

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
         */
        void drawPoint(const Cairo::RefPtr<Cairo::Context> &cr, const Cairo::Matrix &trans,
                double x, double y, const PointType &type);

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

    private:
        // Reference to the parent
        GUI &gui_;
        // simulation object
        std::shared_ptr<eris::Simulation> sim_;
        // The bounds of the graph
        std::array<double, 4> bounds_;
        // Constants for access into bounds_
        static constexpr size_t TOP = 0, RIGHT = 1, BOTTOM = 2, LEFT = 3;
};

}
