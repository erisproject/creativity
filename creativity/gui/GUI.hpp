#pragma once
#include <eris/Simulation.hpp>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <memory>
#include <glibmm/dispatcher.h>
#include <gtkmm.h>
#include <eris/noncopyable.hpp>
#include <unordered_map>
#include <array>
#include <list>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include "creativity/gui/GraphArea.hpp"
#include "creativity/gui/InfoWindow.hpp"
#include "creativity/gui/ReaderCols.hpp"
#include "creativity/gui/ReaderStore.hpp"
#include "creativity/gui/BookStore.hpp"

namespace sigc { SIGC_FUNCTORS_DEDUCE_RESULT_TYPE_WITH_DECLTYPE }

namespace creativity { namespace gui {

/** Class that runs a GUI in a thread, collecting its events into a queue to be processed
 * periodically (e.g. between iterations) via GUIShim.
 */
class GUI : eris::noncopyable {
    public:
        struct Parameter; // forward declaration

        /** Creates a new GUI object.  The GUI is not set up and started until the start() method is
         * called.
         */
        GUI(
                /** The simulation object this GUI is for; the simulation is not expected to be
                 * populated yet. Typically called with an Eris<Simulation> object (which is
                 * castable to an std::shared_ptr) */
                std::shared_ptr<eris::Simulation> eris,
                /** A function to call with the GUI simulation parameters (in a series of
                 * GUI.Parameter structs) when the user configures the simulation via the GUI.  This
                 * function should set up the simulation but not start it: it will be followed
                 * immediately followed by a `run` call (assuming no error occurs).
                 *
                 * If the simulation cannot be started (for example, because some parameters are
                 * invalid) the function should throw a GUI::Exception; the .what() value of the
                 * exception will be displayed to the user as an error message. */
                std::function<void(Parameter param)> setup,
                /** Called to run the simulation for count periods. */
                std::function<void(unsigned int rounds)> run,
                /** Called when the user hits the stop button in the GUI. This should pause the
                 * current simulation run loop until either resume or run are called. */
                std::function<void()> stop,
                /** Called when the user resumes the simulation in the GUI after having stopped it.
                 * This should resume running of the current simulation run loop (if there are
                 * iterations remaining). */
                std::function<void()> resume,
                /** Called when the user wishes to step forward one period then stop. */
                std::function<void()> step,
                /** Called when the user quits the GUI. */
                std::function<void()> quit
        );

        // Destructor: tells the thread to quit and then waits until it does
        virtual ~GUI();

        /** Starts the GUI, reading the glade file and starting the GUI thread.
         */
        void start(int argc, char *argv[]);

        /// Simple wrapper around std::runtime_error used to send error messages back to the user.
        class Exception : public std::runtime_error {
            public:
                /** Constructor: takes an exception message to display to the user. */
                Exception(const std::string &what);
        };

        /** Checks whether the GUI has generated any events and, if so, processes them.  If there
         * are no pending events, this returns immediately.
         */
        void checkEvents();

        /** If the GUI has queued events, processes them.  If there are no queued events, this waits
         * for an event to occur then processes it.
         */
        void waitEvents();

        /** Signals the GUI thread that the visual objects have been updated and the graph needs to
         * be redrawn.  Note that this is not, by default, synchronous: if the GUI thread is
         * currently busy, the redraw might not actually happen until it becomes idle.  It is also
         * collapsed: if multiple redraw() signals are queued, the thread ignores all but the last
         * one.
         *
         * If the optional parameter is given and true, this blocks until the thread finishes a
         * redraw.
         */
        void redraw(bool sync=false);

        /** Sends a signal to the GUI thread that the simulation setup is complete.  The GUI will
         * disable the various simulation setup options and start/initialize buttons and switch to
         * the visualization tab. */
        void initialized();

        /** Sends a signal to the GUI thread that the simulation has started or resumed running. */
        void running();

        /** Signals the GUI thread of a progress update.
         * \param end the last simulation stage of the current run
         * \param speed the speed (in iterations per second)
         * Other progress variables are read from the simulation.
         */
        void progress(const unsigned long &end, const double &speed);

        /** Sends a signal to the GUI thread that the simulation has stopped running.
         *
         * \param manual is true if this stopped because the user hit "Stop" during an active run, false
         * otherwise (e.g. initialization without running, or the configured number of steps completed).
         */
        void stopped(bool manual);

        /** Sends a signal to the GUI thread that an error has occurred.  Takes a string to display
         * in a dialog to the user. */
        void error(std::string message);

        /// Parameter types for a Parameter.
        enum class ParamType {
            /// Sets the number of dimensions in `.ul`
            dimensions,
            /// Sets the number of readers in `.ul`
            readers,
            /// Sets the book location standard deviation in `.dbl`
            book_sd,
            /// Sets the perceived quality standard deviation in `.dbl`
            quality_draw_sd,
            /// Sets the fixed cost of keeping a book on the market in `.dbl`
            cost_fixed,
            /// Sets the per-unit cost of copies of a book in `.dbl`
            cost_unit,
            /// Sets the simulation speed limit in `.dur_ms`
            speed_limit,
            /// Sets the minimum redraw period in `.dur_ms`
            redraw,
            /// Number of threads to use in `.ul`
            threads,
            /// Sent by the GUI to indicate that some parameters are being changed.
            begin,
            /// Fired when setup ends unsuccessfully because when one or more setup parameters threw exceptions
            erred,
            /// Fired when setup ends successfully (no setup parameter threw an exception)
            finished
        };
        /** The parameter struct for passing a configured value back from the GUI.  The GUI always
         * sends a `begin` followed by zero or more settings then either `erred` or `finished` (the
         * former if one or more exceptions occurred in the settings, the latter if no exceptions
         * occurred).
         */
        struct Parameter {
            /// The parameter type
            ParamType param;
            /// A value of various types
            union {
                bool bl;
                unsigned long ul;
                long l;
                unsigned int ui;
                int i;
                double dbl;
                std::chrono::milliseconds dur_ms;
            };
        };

        /** Simple class containing event information.
         */
        class Event {
            public:
                /// The types of events that the GUI thread can send to the main thread
                enum class Type {
                    setup, ///< An event triggered when the user hits "begin" with configuration data
                    run, ///< An event triggered when starting a new or resuming a stopped simulation
                    stop, ///< The user hit the "stop" button to pause the simulation.
                    resume, ///< The user hit the "resume" button to unpause the simulation.
                    step, ///< The user hit the "step" button to unpause the simulation for one step.
                    quit, ///< Sent when the user quits the application
                    redraw_complete ///< Sent to indicate that a redraw is complete.  Used internally for redraw(true) calls.
                };

                /// Will be true if this is a fake Event indicating that no events are pending.
                bool none;

                /** True if this is actually an event, false if this reflects that there are no
                 * events pending.
                 */
                operator bool();

                /** The type of the event.  Only valid if the Event object evaluates in boolean
                 * context to true (or, equivalently, none is false).
                 */
                Type type;

                /** Will be set to a list of configured GUI::Parameters if this is a setup event.
                 * The `begin` and `finished` GUI::Parameter meta-values are *not* included.
                 */
                std::vector<Parameter> parameters;

                /** unsigned long integer value associated with the event.  Used by Event::Type::run
                 * to send the number of iterations to run.
                 */
                unsigned long ul;

            private:
                friend class GUI;
                /// Default constructor: the "no events" event
                Event();
                /// Constructs an event with no parameters or value
                Event(Type t);
                /// Constructs an event with Parameters.  The vector should not contain
                /// begin/finished parameter events; these will be sent appropriately.
                Event(Type t, std::vector<Parameter> &&p);
                /// Constructs an event with an unsigned long value
                Event(Type t, const unsigned long &ul);
        };

        /** Processes events sent by the GUI until at least one event of the given type has been
         * received.
         */
        void waitForEvent(Event::Type t);

        /** Utility method used by various GUI classes to convert a Position to a user-displayable string.
         */
        static std::string pos_to_string(const eris::Position &pos);

    private:
        /// The thread the GUI is running in.  Set during construction.
        std::thread gui_thread_{};

        /** The GTK::Application.  Created in start(); should only be touched after startup by the
         * GUI thread.
         */
        decltype(Gtk::Application::create()) app_;

        /** The Gtk::Builder that creates the GUI from the glade file.  Should only be touched after
         * startup by the GUI thread.
         */
        decltype(Gtk::Builder::create()) builder_;

        /** The main window. */
        std::shared_ptr<Gtk::Window> main_window_;

        /// Will be true once the thread has finished setting itself up and started its mainloop.
        bool thread_running_{false};

        /// Command-line parameters can override the gui.glade default values for settings
        std::unordered_map<std::string, double> default_override_;

        /** Dispatcher for receiving a signal inside the GUI thread.  When this gets sent, the GUI
         * checks signal_queue_ for instructions and acts accordingly.
         */
        std::unique_ptr<Glib::Dispatcher> dispatcher_;

        /** Mutex controlling access to the data shared between the main and GUI threads.  This
         * mutex also handles synchronization during thread startup and guards access to the queue
         * of events (setup, run, etc.) emitted by the GUI and status variable accessed by the
         * GUI.  It also guards dispatcher_, which the thread deletes when it quits.
         */
        std::mutex mutex_;

        /// Used to signal the main thread that the initialization is done
        std::condition_variable cv_;

        /** Swaps out the queue elements into a local queue, unlocks the lock, then processes the
         * swapped out queue elements.  This is the common code for waitEvents() and checkEvents().
         */
        void processEvents_(std::unique_lock<decltype(mutex_)> &lock);

        /** Signal class to send a signal with a tuple of arbitrary data. */
        class Signal {
            public:
                /// The types of signal that can be sent to the GUI thread
                enum class Type { redraw, initialized, running, progress, stopped, error, quit };

                /// The type of signal
                Type type;

                /// Possible string associated with the signal.
                std::string message;
                /// Possible bool
                bool boolean;
                /// Possible unsigned longs associated with the signal.
                std::vector<unsigned long> uls;
                /// Doubles
                std::vector<double> doubles;

                Signal() = delete;
                /// Signal without data (NB: this implicitly means any Signal::Type can be cast to a Signal):
                Signal(Type type) : type(type) {}
                /// Signal with a message:
                Signal(Type type, const std::string &message) : type{type}, message{message} {}
                /// Signal with boolean
                Signal(Type type, bool b) : type{type}, boolean{b} {}
        };

        /// The queue of signals that have been sent to but not yet processed by the GUI thread
        std::list<Signal> signal_queue_;

        /** Queues a signal for the GUI thread and signals it. */
        void queueSignal(Signal &&s);

        /// The queue of events that have been sent by the GUI thread but not yet processed
        std::list<Event> event_queue_;

        /** Constructs and adds an Event to the event_queue_, taking care of mutex locking and cv_
         * signalling.
         */
        template <typename... Args>
        void queueEvent(Args&&... args);

        friend class GraphArea;

        /// Called by the GUI thread to process the signal queue
        void thr_signal();

        /** Started in a thread; sets up the various graphical stuff and starts the main loop.
         * `app_` and `builder_` must be set up in the main thread (so that builder errors happen
         * *before* starting the thread).
         */
        void thr_run();

        /** Updates the reader information on the "Agents" tab.
         */
        void thr_update_readers();

        /** Updates the book information on the "Books" tab.
         */
        void thr_update_books();

        /** Opens a dialog for the given member, which must be either a Reader or a Book.  If the
         * dialog is already open, it is presented again (which is window manager dependent, but
         * generally means bringing to the top and/or focussing). */
        void thr_info_dialog(eris::SharedMember<eris::Member> member);

        /** Sets up the simulation based on the current GUI parameter values.
         */
        void setupSim();

        /** Starts the simulation (for however many periods specified in the "set_periods" GUI
         * configuration option).  Must be preceeded by a call to setupSim();
         */
        void runSim();

        /** The custom graph area */
        std::unique_ptr<GraphArea> graph_;

        /** The various objects used for the Agents tab */
        std::unique_ptr<ReaderCols> rdr_cols_;
        Glib::RefPtr<ReaderStore> rdr_list_;
        std::unique_ptr<Gtk::TreeView> rdr_tree_;
        // The latest reader added (since new eris_id_t values are strictly increasing, we can use
        // this to identify new readers).
        eris::eris_id_t rdr_biggest_id_ = 0;

        /** The various objects used for the Books tab */
        Glib::RefPtr<BookStore> bk_model_;
        std::unique_ptr<Gtk::TreeView> bk_tree_;

        /** Obtains a widget from the current Gtk::Builder and returns it (as a pointer).
         */
        template <class T>
        T* widget(const std::string &widget_name);

        /** Returns the value of a Gtk::SpinButton. */
        double sb(const std::string &widget_name);

        /** Returns the value of a Gtk::SpinButton as an int. */
        int sb_int(const std::string &widget_name);

        /// The simulation object
        std::shared_ptr<eris::Simulation> sim_;
        // The callbacks for GUI thread events
        std::function<void(GUI::Parameter)> on_setup_;
        std::function<void(unsigned int count)> on_run_;
        std::function<void()> on_stop_;
        std::function<void()> on_resume_;
        std::function<void()> on_step_;
        std::function<void()> on_quit_;

        /** Handles a single event received from the GUI thread. */
        void handleEvent(const GUI::Event &event);

        /** Currently open reader/book dialogs */
        std::unordered_map<eris::eris_id_t, InfoWindow> info_windows_;

        typedef boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian> rt_point;
        typedef std::pair<rt_point, eris::SharedMember<eris::Member>> rt_val;
        /** An rtree (regenerated each time the simulation stops/pauses) of reader/book points. */
        boost::geometry::index::rtree<rt_val, boost::geometry::index::rstar<16>> rtree_;

        std::function<bool(GdkEventMotion* event)> motion_handler_;
        sigc::connection motion_handler_conn_;
        Glib::RefPtr<Gdk::Cursor> hand_;
};

template <typename... Args>
void GUI::queueEvent(Args&&... args) {
    std::unique_lock<std::mutex> lock(mutex_);
    event_queue_.push_back(Event(std::forward<Args>(args)...));
    lock.unlock();
    cv_.notify_all();
}

} }
