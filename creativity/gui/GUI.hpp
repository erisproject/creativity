#pragma once
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
#include "creativity/state/State.hpp"
#include "creativity/state/Storage.hpp"

namespace sigc { SIGC_FUNCTORS_DEDUCE_RESULT_TYPE_WITH_DECLTYPE }

namespace creativity { namespace gui {

// Forward declarations
class ReaderStore;
class BookStore;

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
                /** The Creativity object containing the simulation storage.  The object does not
                 * need to be initialized (by calling setup()): in particular, if the GUI is
                 * displaying a previously stored simulation run, it typically won't be.
                 */
                std::shared_ptr<Creativity> creativity,

                /** A function to call with the GUI simulation parameters (in a series of
                 * GUI.Parameter structs) when the user configures the simulation via the GUI.  This
                 * function should set up the simulation but not start it: it will be followed
                 * immediately followed by a `run` call (assuming no error occurs).
                 *
                 * If the simulation cannot be started (for example, because some parameters are
                 * invalid) the function should throw a GUI::Exception; the .what() value of the
                 * exception will be displayed to the user as an error message. */
                std::function<void(Parameter param)> setup,

                /** Called to run the simulation for `rounds` periods. */
                std::function<void(eris::eris_time_t rounds)> run,

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

        /** Signals the GUI thread that one or more states have been completed.  The GUI will record
         * the new number of states.  If the user was viewing the most recent state, the GUI will be
         * updated to display the new most recent state, otherwise the GUI remains on the current
         * state.
         *
         * The optional parameter allows the main thread to instruct the GUI to switch to the given
         * state; if omitted, the GUI switches to the last state if the user was already there
         * before the signal, and otherwise doesn't change the state.
         */
        void newStates(unsigned long switch_to = (unsigned long) -1);

        /** Sends a signal to the GUI thread that the simulation setup is complete.  The GUI will
         * disable the various simulation setup options and start/initialize buttons and switch to
         * the visualization tab.  The `states` variable given to the constructor must have been
         * initialized with the initial simulation state: this method sends a new_states signal, to
         * let the GUI thread know that there is a new state. */
        void initialized();

        /** Sends a signal to the GUI thread that the simulation has started or resumed running. */
        void running();

        /** Signals the GUI thread of a progress update.
         * \param t the current simulation state
         * \param end the simulation stage at which the current run will stop
         * \param speed the speed (in iterations per second)
         */
        void progress(unsigned long t, unsigned long end, double speed);

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
            save_as, ///< The file to save simulation data to (will be overwritten)
            load, ///< The file to load existing simulation data from
            seed, ///< Sets the seed value for eris::Random::seed in `.ul`
            threads, ///< Number of threads to use in `.ul`
            begin, ///< Sent by the GUI to indicate that some parameters are being changed.
            erred, ///< Fired when setup ends unsuccessfully because when one or more setup parameters threw exceptions
            finished ///< Fired when setup ends successfully (no setup parameter threw an exception)
        };
        /** The parameter struct for passing a configured value back from the GUI.  The GUI always
         * sends a `begin` followed by zero or more settings then either `erred` or `finished` (the
         * former if one or more exceptions occurred in the settings, the latter if no exceptions
         * occurred).
         */
        struct Parameter {
            /// The parameter type
            ParamType param;
            /// A value of various primitive types
            union {
                bool bl;
                unsigned long ul;
                long l;
                unsigned int ui;
                int i;
                double dbl;
                std::chrono::milliseconds dur_ms;
                void *ptr;
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
                    none ///< Fake Event indicating that no events are pending
                };

                /** True if this is actually an event, false if this reflects that there are no
                 * events pending (that is, an event with `type = Event::Type::none`).
                 */
                operator bool();

                /// The type of the event.
                Type type = Type::none;

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
                Event() = default;
                /// Constructs an event with no parameters or value
                Event(Type t);
                /// Constructs an event with Parameters.  The vector should not contain
                /// begin/finished parameter events; these will be sent appropriately.
                Event(Type t, std::vector<Parameter> &&p);
                /// Constructs an event with an unsigned long value
                Event(Type t, unsigned long ul);
        };

        /** Processes events sent by the GUI until at least one event of the given type has been
         * received.
         */
        void waitForEvent(Event::Type t);

        /** Utility method used by various GUI classes to convert a Position to a user-displayable string.
         */
        static std::string pos_to_string(const eris::Position &pos);

    private:
        /// The Creativity object
        std::shared_ptr<Creativity> creativity_;

        /// The thread the GUI is running in.  Set during construction.
        std::thread gui_thread_{};

        /** The GTK::Application.  Created in start(); should only be touched after startup by the
         * GUI thread.
         */
        Glib::RefPtr<Gtk::Application> app_;

        /** The Gtk::Builder that creates the GUI from the glade file.  Should only be touched after
         * startup by the GUI thread.
         */
        Glib::RefPtr<Gtk::Builder> builder_;

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

        /** Signal class to send a signal with various types of data from the main thread to the GUI. */
        class Signal {
            public:
                /// The types of signal that can be sent to the GUI thread
                enum class Type {
                    /** The main thread may have added new states to the state storage; the GUI
                     * should lock it and update associated internal variables if necessary.  This
                     * state may optionally includes a single value in `.uls` indicating a state
                     * number the GUI should switch to; if this value is omitted, the GUI will
                     * switch to the last state if the user was already on the last state, and
                     * otherwise not change the current state.
                     */
                    new_states,
                    initialized, ///< The main thread is finished with initial setup
                    running, ///< The main thread has started or resumed running
                    progress, ///< Updates the current progress in the GUI
                    stopped, ///< The main thread has stopped
                    error, ///< The main thread encountered an error (typically during setup)
                    quit, ///< The GUI object is being destroyed and should terminate
                    none ///< Placeholder signal, set when a Signal is default constructed
                };

                /// The type of signal
                Type type = Type::none;

                /// Possible string associated with the signal.
                std::string message;
                /// Possible bool
                bool boolean;
                /// Possible unsigned longs associated with the signal.
                std::vector<unsigned long> uls;
                /// Doubles
                std::vector<double> doubles;

                /// Default constructor: constructs a signal with type set to Signal::Type::none
                Signal() = default;
                /// Signal without data (NB: this implicitly means any Signal::Type can be cast to a Signal):
                Signal(Type type) : type(type) {}
                /// Signal with a message:
                Signal(Type type, const std::string &message) : type{type}, message{message} {}
                /// Signal with boolean
                Signal(Type type, bool b) : type{type}, boolean{b} {}
                /// Signal with a single double value
                Signal(Type type, double d) : type{type}, doubles(1, d) {}
                /// Signal with a single unsigned long value
                Signal(Type type, unsigned long ul) : type{type}, uls(1, ul) {}
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

        /** When called, this updates the simulation parameters displayed in the GUI to match the
         * current creativity_->parameters values.
         */
        void thr_update_parameters();

        /** Called to update the GUI tabs to display the state at the given time.
         */
        void thr_set_state(unsigned long t);

        /** The state currently be displayed. */
        unsigned long state_curr_ = (unsigned long) -1;

        /** The current number of states known to the GUI.  `creativity_.storage.size() >=
         * state_num_` is guaranteeded to be true.  There may be more in the `creativity_.storage`
         * variable if the GUI hasn't received and processed the notification yet.
         */
        unsigned long state_num_ = 0;

        /** Opens a dialog for the given member, which must be either a Reader or a Book.  If the
         * dialog is already open, it is presented again (which is window manager dependent, but
         * generally means bringing to the top and/or focussing).
         *
         * \throws std::out_of_range if `member_id` does not exist in the current state
         */
        void thr_info_dialog(eris::eris_id_t member_id);

        /** Loads a previously-run simulation from the currently selected file.
         */
        void loadSim();

        /** Sets up a new simulation based on the current GUI parameter values.
         */
        void setupSim();

        /** Starts the simulation (for however many periods specified in the "set_periods" GUI
         * configuration option).  Must be preceeded by a call to setupSim();
         */
        void runSim();

        /** The custom graph area */
        std::unique_ptr<GraphArea> graph_;

        /** The various objects used for the Agents tab */
        Gtk::ScrolledWindow *rdr_win_;
        std::vector<Glib::RefPtr<ReaderStore>> rdr_models_;
        std::vector<std::unique_ptr<Gtk::TreeView>> rdr_trees_;

        /** The various objects used for the Books tab */
        Gtk::ScrolledWindow *bk_win_;
        std::vector<Glib::RefPtr<BookStore>> bk_models_;
        std::vector<std::unique_ptr<Gtk::TreeView>> bk_trees_;

        /** Obtains a widget from the current Gtk::Builder and returns it (as a pointer).
         */
        template <class T>
        typename std::enable_if<std::is_base_of<Gtk::Widget, T>::value, T*>::type
        widget(const std::string &widget_name) {
            T* widget;
            builder_->get_widget(widget_name, widget);
            return widget;
        }

        /** Returns the value of a Gtk::SpinButton. */
        double sb(const std::string &widget_name);

        /** Returns the value of a Gtk::SpinButton as an int. */
        int sb_int(const std::string &widget_name);

        /** Storage for the "save as" and "load" filenames; pointers to these are sent back when
         * needed to the main thread.
         */
        std::string save_, load_;

        // The file filter.  Access via fileFilter().
        mutable Glib::RefPtr<Gtk::FileFilter> ff_;

        /** Returns a file filter for *.crstate files. */
        decltype(ff_) fileFilter() const;

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
        typedef std::pair<rt_point, eris::eris_id_t> rt_val;
        using RTree = boost::geometry::index::rtree<rt_val, boost::geometry::index::rstar<16>>;
        /** rtrees of reader/book points for each state. */
        std::vector<std::unique_ptr<RTree>> rtrees_;

        /** Returns the `n` nearest points/eris_id_t pairs to the given rt_point.  May returns fewer
         * than `n` (including 0) if there are not `n` points within `radius` pixels.
         *
         * \param point the graph point (not canvas point) around which to search
         * \param n the maximum number of rt_vals to return; the first element is the closest point.
         * Defaults to 1.
         */
        std::vector<rt_val> thr_nearest(const rt_point &point, int n = 1);

        /** Called when the visibility of readers, live books, or dead books have changed to throw
         * out all cached rtrees and rebuild the current state's rtree.
         */
        void thr_reset_rtrees();

        /** Given a newly-constructed, empty rtree and a state, this populates the rtree with the
         * visible books and readers of the state.
         */
        void thr_init_rtree(RTree &rt, const std::shared_ptr<const state::State> &state) const;

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
