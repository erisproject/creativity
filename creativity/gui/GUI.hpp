#pragma once
#include "creativity/cmdargs/GUI.hpp"
#include "creativity/gui/GraphArea.hpp"
#include <eris/noncopyable.hpp>
#include <eris/Position.hpp>
#include <eris/types.hpp>
#include <glibmm/refptr.h>
#include <sigc++/connection.h>
#include <gtkmm/widget.h>
#include <gtkmm/builder.h> // IWYU pragma: keep
#include <boost/geometry/core/cs.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <unordered_map>
#include <list>
#include <chrono>
#include <functional>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <memory>

namespace Gdk { class Cursor; }
namespace Glib { class Dispatcher; }
namespace Gtk { class Application; }
namespace Gtk { class FileFilter; }
namespace Gtk { class ScrolledWindow; }
namespace Gtk { class TreeView; }
namespace Gtk { class Window; }
namespace creativity { class Creativity; }
namespace creativity { namespace state { class State; } }

namespace creativity {
/// Namespace for all gtkmm creativity GUI classes
namespace gui {

// Forward declarations
class ReaderStore;
class BookStore;
class InfoWindow;

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
                 *
                 * Note that the creativity object must persist as long as the created GUI object.
                 */
                Creativity &creativity,

                /** A function to call with GUI simulation parameters when the user configures the
                 * simulation via the GUI.  Note that this doesn't include simulation settings (like
                 * the number of readers)--those are set directly in the settings object.
                 *
                 * If the simulation cannot be started (for example, because some parameters are
                 * invalid) the function should throw an exception derived from std::exception; the
                 * .what() value of the exception will be displayed to the user as an error message.
                 */
                std::function<void(Parameter param)> configure,

                /** A function to call to initialize the simulation without starting it.  If
                 * something goes wrong, this should throw an exception derived from std::exception:
                 * the `.what()` value of the exception will be displayed to the user as an error
                 * message.
                 */
                std::function<void()> initialize,

                /** A function called to specify the number of simulation periods, and called again
                 * whenever the user changes the number of simulation periods. */
                std::function<void(eris::eris_time_t end)> change_periods,

                /** Called to run or continue running the simulation until the number of periods
                 * sent by the last `change_periods` call have been completed.  `initialize` will be
                 * called before the first call to this, and `change_periods` will have been called
                 * one or more times.
                 */
                std::function<void()> run,

                /** Called when the user hits the stop button in the GUI. This should pause the
                 * current simulation run loop until run or step are called. */
                std::function<void()> stop,

                /** Called when the user wishes to step forward one period then stop. */
                std::function<void()> step,

                /** Called when the user quits the GUI. */
                std::function<void()> quit
        );

        // Destructor: tells the thread to quit and then waits until it does
        virtual ~GUI();

        /** Starts the GUI, reading the glade file and starting the GUI thread.
         */
        void start(const cmdargs::GUI &args);

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
        void newStates(eris::eris_time_t switch_to = (eris::eris_time_t) -1);

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
        void progress(eris::eris_time_t t, eris::eris_time_t end, double speed);

        /** Sends a signal to the GUI thread that the simulation has stopped running.
         *
         * \param more is true if there are still more simulation periods to run (as specified in
         * the Periods setting).
         */
        void stopped(bool more);

        /** Sends a signal to the GUI thread that an error has occurred.  Takes a string to display
         * in a dialog to the user. */
        void error(std::string message);

        /// Parameter types for a Parameter.
        enum class ParamType {
            save_as, ///< The file to save simulation data to (will be overwritten)
            load, ///< The file to load existing simulation data from
            seed, ///< Sets the seed value for eris::Random::seed in `.ul`
            threads ///< Number of threads to use in `.ul`
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
                    configure, ///< An event that configures a Parameter
                    initialize, ///< An event triggered when the user hits "begin" with configuration data
                    periods, ///< An event indicating that the number of periods is now known or has changed
                    run, ///< An event triggered when starting a new or resuming a stopped simulation
                    stop, ///< The user hit the "stop" button to pause the simulation.
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
                Parameter parameter;

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
                /// Constructs a parameter event with the given Parameter.
                Event(const Parameter &p);
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
        Creativity &creativity_;

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

        /** The main window title; this comes from the glade file, but appends the version. */
        std::string main_window_title_;

        /// Will be true once the thread has finished setting itself up and started its mainloop.
        bool thread_running_{false};

        /// Whether the piracy tick has been added to the period slider
        bool piracy_tick_added_ = false;

        /// Whether the public tick has been added to the period slider
        bool public_tick_added_ = false;

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
                /// Possible uint64_t values associated with the signal.
                std::vector<uint64_t> uls;
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
                /// Signal with a single unsigned integer value (value will be stored in uls)
                template<typename UIntT, typename = typename std::enable_if<std::is_integral<UIntT>::value and std::is_unsigned<UIntT>::value>::type>
                Signal(Type type, UIntT ul) : type{type}, uls(1, ul) {}
                /// Signal with a single uint32_t value (will be stored in uls)
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
         *
         * The given cmdargs::GUI object should have already parsed command-line arguments.
         */
        void thr_run(const cmdargs::GUI &args);

        /** Sets up the handlers for visualization setting `field`, which is controlled by
         * "enable_<field>" and "colour_<field>" GUI elements.  To set up just one or the other of
         * these, see thr_connect_vis_colour() or thr_connect_vis_enabled() -- this method simply
         * calls both of those.
         *
         * \param the name of the field
         * \param enabled the bool reference controlling whether this field is enabled
         * \param colour the Colour reference controlling the colour associated with this field
         * \param affects_rtree whether changing this setting requires a rebuild of the rtree
         * (mapping clickable elements): this should be true if the field enables/disables clickable
         * elements.
         * \param needs_all if non-empty, this field will be set insensitive if any of the given
         * elements are disabled.
         * \param needs_any if non-empty, this field will be set insensitive if all of the given
         * elements are disabled.
         */
        void thr_connect_vis_setting(
                const std::string &field,
                bool &enabled,
                GraphArea::Colour &colour,
                bool affects_rtree = false,
                const std::vector<std::string> &needs_all = {},
                const std::vector<std::string> &needs_any = {});

        /** Sets up the "enable_field" handler in the GUI to adjust the internal GUI parameter
         * enabling or disabling the feature.
         *
         * \param the field name (not prefixed with "enable_")
         * \param enabled reference to the bool value to enable/disable as the GUI setting is
         * changed
         * \param affects_rtree whether changing this setting requires a rebuild of the rtree
         * (mapping clickable elements): this should be true if the field enables/disables clickable
         * elements.
         * \param needs_all if non-empty, this field will be set insensitive if any of the given
         * elements are disabled.
         * \param needs_any if non-empty, this field will be set insensitive if all of the given
         * elements are disabled.
         * \sa thr_connect_vis_setting
         */
        void thr_connect_vis_enabled(
                const std::string &field,
                bool &enabled,
                bool affects_rtree = false,
                const std::vector<std::string> &needs_all = {},
                const std::vector<std::string> &needs_any = {});

        /** Sets up the "colour_field" handler in the GUI to adjust the internal GUI colour for the
         * associated field.
         *
         * \param the field name (not prefixed with "enable_")
         * \param colour the internal Colour value to change when the GUI setting is changed.
         * \param needs_all if non-empty, this field will be set insensitive if any of the given
         * elements are disabled.
         * \param needs_any if non-empty, this field will be set insensitive if all of the given
         * elements are disabled.
         * \sa thr_connect_vis_setting
         */
        void thr_connect_vis_colour(
                const std::string &field,
                GraphArea::Colour &colour,
                const std::vector<std::string> &needs_all = {},
                const std::vector<std::string> &needs_any = {});

        /** Sets the widget `widget` to be enabled whenever everything in needs_all is enabled, and
         * anything in needs_any is enabled.  (Both are considered satisfied when the relevent
         * vector is empty).  Does nothing at all if both are empty.
         */
        void thr_connect_vis_deps(
                Gtk::Widget *widget,
                const std::vector<std::string> &needs_all,
                const std::vector<std::string> &needs_any);

        /** When called, this updates the simulation parameters displayed in the GUI to match the
         * current creativity_.parameters values.
         */
        void thr_update_parameters();

        /** Called to update the GUI tabs to display the state at the given time.
         */
        void thr_set_state(eris::eris_time_t t);

        /** The state currently be displayed. */
        eris::eris_time_t state_curr_ = (eris::eris_time_t) -1;

        /** The current number of states known to the GUI.  `creativity_.storage.size() >=
         * state_num_` is guaranteeded to be true.  There may be more in the `creativity_.storage`
         * variable if the GUI hasn't received and processed the notification yet.
         */
        eris::eris_time_t state_num_ = 0;

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

        /** Writes the current simulation and any future simulation states to the selected file.
         */
        void saveSim();

        /** Sets up a new simulation based on the current GUI parameter values.
         */
        void initializeSim();

        /** Starts the simulation (for however many periods specified in the "set_periods" GUI
         * configuration option).  Must be preceeded by a call to initializeSim();
         */
        void runSim();

        /** The custom graph area */
        std::unique_ptr<GraphArea> graph_;

        /** The Agents tab */
        Gtk::ScrolledWindow *rdr_win_;

        /** The Books tab */
        Gtk::ScrolledWindow *bk_win_;

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

        /** Returns a file filter for *.crstate and *.crstate.xz files. */
        decltype(ff_) fileFilter() const;

        // The callbacks for GUI thread events
        std::function<void(GUI::Parameter)> on_configure_;
        std::function<void()> on_initialize_;
        std::function<void(unsigned int count)> on_change_periods_;
        std::function<void()> on_run_;
        std::function<void()> on_stop_;
        std::function<void()> on_step_;
        std::function<void()> on_quit_;

        /** Handles a single event received from the GUI thread. */
        void handleEvent(const GUI::Event &event);

        /** Currently open reader/book dialogs */
        std::unordered_map<eris::eris_id_t, std::unique_ptr<InfoWindow>> info_windows_;

        using rt_point = boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian>;
        using rt_val = std::pair<rt_point, eris::eris_id_t>;
        using RTree = boost::geometry::index::rtree<rt_val, boost::geometry::index::rstar<16>>;

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

        /** The various per-period data objects, used for constructing the data elements for a
         * period */
        struct temporal_data {
            Glib::RefPtr<ReaderStore> rdr_model;
            std::unique_ptr<Gtk::TreeView> rdr_tree;
            Glib::RefPtr<BookStore> bk_model;
            std::unique_ptr<Gtk::TreeView> bk_tree;
            std::unique_ptr<RTree> rtree;
        };

        /// Cache of most-recently-used temporal data
        std::list<std::pair<eris::eris_time_t, temporal_data>> temporal_cache_;
        /// Size of the temporal cache data; also used by GraphArea for the size of its temporal
        /// image data cache
        unsigned int temporal_cache_size_ = 10;

};

template <typename... Args>
void GUI::queueEvent(Args&&... args) {
    std::unique_lock<std::mutex> lock(mutex_);
    event_queue_.push_back(Event(std::forward<Args>(args)...));
    lock.unlock();
    cv_.notify_all();
}

} }
