#pragma once
#include <thread>
#include <mutex>
#include <condition_variable>
#include <memory>
#include <sigc++/sigc++.h>
#include <glibmm/dispatcher.h>
#include <gtkmm.h>
#include <boost/utility.hpp>
#include <unordered_map>
#include <array>

namespace sigc { SIGC_FUNCTORS_DEDUCE_RESULT_TYPE_WITH_DECLTYPE }

namespace creativity {

/** Class that runs a GUI in a thread, collecting its events into a queue to be processed
 * periodically (e.g. between iterations) via GUIShim.
 */
class GUI : boost::noncopyable {
    public:
        /// The parameters configurable through the GUI.
        struct Parameters {
            unsigned int dimensions;
            unsigned int readers;
            unsigned int publishers;
        };

        /// Simple wrapper around std::runtime_error used to send error messages back to the user.
        class Exception : public std::runtime_error {
            public:
                Exception(const std::string &what);
        };

        /** Creates a new GUI object.  The GUI is not set up and started until the start() method is
         * called.
         */
        GUI();

        /** Starts the GUI, reading the glade file and starting the GUI thread.
         */
        void start(int argc, char *argv[]);

        /** Sends a signal to the GUI thread to quit, and waits for the thread to exit before
         * returning.
         */
        void quit();

        /** Point types supported by addPoint */
        enum class PointType {
            X, //!< A diagonal cross
            CROSS, //!< A horizontal/vertical cross
            SQUARE //!< A small square
        };
        /** Adds a point to the graph at the given position. Returns a value that can later be
         * passed to removePoint to clear the point. Does *not* update the GUI visualization with
         * the new point: for that, call redraw(). */
        unsigned long addPoint(const double &x, const double &y, const PointType &type);

        /** Removes a point.  Takes a point id as returned by addPoint(). Returns true if the point
         * was removed, false if the point wasn't on the graph. Does *not* update the GUI until
         * redraw() is called. */
        bool removePoint(const unsigned long &id);

        /** Clears all points on the graph. Does *not* update the GUI until redraw() is called. */
        void clearPoints();

        /** Circle types supported by addCircle() */
        enum class CircleType {
            A, //!< Circle style A
            B //!< Circle style B, will look different from A
        };
        /** Adds a circle to the graph at the given center and radius. Returns an id that can later
         * be passed to removeCircle to remove it. Does *not* update the GUI until redraw() is
         * called. */
        unsigned long addCircle(const double &cx, const double &cy, const double &r, const CircleType &type);

        /** Removes a circle.  Takes a circle id as returned by addCircle(). Returns true if the
         * circle was removed, false if the circle wasn't on the graph. Does *not* update the GUI
         * until redraw() is called. */
        bool removeCircle(const unsigned long &id);

        /** Clears all circles on the graph. Does *not* update the GUI until redraw() is called. */
        void clearCircles();

        /** The radius of markers representing points, in pixels.  For example, a value of 5 means
         * that CROSS points will have horizontal and vertical lines extending a distance of 5,
         * while X points will have diagonal lines with a length of 5.  For SQUARE points the
         * diagonals of the square will have length 10 (and so edges will be \f$\5 \sqrt{2}\f$
         * pixels long).
         */
        double pointMarkerRadius = 5.0;

        /** Signals the GUI thread that the visual objects have been updated and the graph needs to
         * be redrawn.  Note that this is not synchronous: if the GUI thread is currently busy, the
         * redraw might not actually happen until it becomes idle.
         */
        void redraw();

    protected:
        /// The thread the GUI is running in.  Set during construction.
        std::unique_ptr<std::thread> gui_thread_;

        /** Dispatcher for receiving a signal inside the GUI thread that the visual display needs to
         * be updated.  Set during thread creation.  Note that this is thread-safe: you do not need
         * to obtain a mutex lock before signaling with it.
         */
        std::unique_ptr<Glib::Dispatcher> redraw_dispatcher_;

        /// Dispatcher used to tell the GUI thread to quit.
        std::unique_ptr<Glib::Dispatcher> quit_dispatcher_;

        /// Will be true once the thread has finished setting itself up and started its mainloop.
        bool thread_running_{false};

        /** Mutex and lock controlling access to the data shared between the main and GUI threads.
         * This lock also handles synchronization during thread startup and guards access to the
         * queue of events (setup, run, thread) emitted by the GUI.
         */
        std::mutex mutex_;

        /// Used to signal the main thread that the initialization is done
        std::condition_variable cv_;

        /** Started in a thread; sets up the various graphical stuff and starts the main loop.
         * Takes the pre-loaded Gtk::Builder (so that builder errors happen *before* starting the
         * thread).
         */
        void thr_run(decltype(Gtk::Application::create()) app, decltype(Gtk::Builder::create()) builder);

        struct Point_ { double x, y; PointType type; };
        struct Circle_ { double cx, cy, r; CircleType type; };
        typedef std::unordered_map<unsigned long, const Point_> pointmap_t;
        typedef std::unordered_map<unsigned long, const Circle_> circlemap_t;
        unsigned long point_id_next_ = 0, circle_id_next_ = 0;
        pointmap_t points_;
        circlemap_t circles_;

        class GraphArea : public Gtk::DrawingArea, boost::noncopyable {
            public:
                /** Creates a graph area that draws in the rectangle bounded by [`bottom', `top'] on the
                 * vertical plane and [`left', `right'] on the horizontal plane.  Specifying a value of
                 * bottom or left larger than top or bottom will flip the respective axis.
                 */
                GraphArea(const double &top, const double &right, const double &bottom, const double &left, GUI &gui);

                /** Draws the current set of points/circles. */
                virtual bool on_draw(const Cairo::RefPtr<Cairo::Context> &cr) override;

                /** Draws a point.  Points have a fixed size that do not depend on the drawing area
                 * size.
                 */
                void drawPoint(const Cairo::RefPtr<Cairo::Context> &cr, const Point_ &point, const Cairo::Matrix &trans);

                /** Draws a circle.  Circles, however, are in graph space, not drawing area space, so
                 * this is actually going to end up drawing ovals (unless the drawing area happens
                 * to be perfectly square).
                 */
                void drawCircle(const Cairo::RefPtr<Cairo::Context> &cr, const Circle_ &circle, const Cairo::Matrix &trans);

            private:
                // Reference to the parent
                GUI &gui_;
                // The bounds of the graph
                const std::array<double, 4> bounds_;
                // Constants for access into bounds_
                static constexpr size_t TOP = 0, RIGHT = 1, BOTTOM = 2, LEFT = 3;
        };

        std::unique_ptr<GraphArea> graph_;

};

}
