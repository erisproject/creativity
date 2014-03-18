#pragma once
#include <eris/Eris.hpp>
#include <eris/Simulation.hpp>
#include <mutex>
#include <memory>
#include <sigc++/sigc++.h>
#include <boost/utility.hpp>
#include "creativity/GUI.hpp"

using namespace eris;

namespace creativity {

/** Class that operates as a shim between the creativity simulation code and the GUI, which handles
 * all graphics-related code in a thread.  This class GUI class maintains its own thread controlling the
 * actual GTK GUI, providing methods to send and receive signal to and from the graphical interface.
 */
class GUIShim : boost::noncopyable {
    public:
        /** Creates a new GUIShim object, which internally creates a new GUI object (starting the
         * GUI in a new thread).
         *
         * Immediately after construction you should call, at least, setupHandler() an
         */
        GUIShim(
                /** The simulation object this GUI is for; the simulation is not expected to be
                 * populated yet. Typically called with an Eris<Simulation> object (which is
                 * castable to an std::shared_ptr) */
                std::shared_ptr<Simulation> eris,
                /** A function to call with the GUI simulation parameters (in a GUI.Parameters
                 * struct) when the user instructs (via the GUI) to start the simulation.  This
                 * function should set up the simulation but not start it: it will be followed
                 * immediately followed by a `resume` call (assuming no error occurs).
                 *
                 * If the simulation cannot be started (for example, because some parameters are
                 * invalid) the function should throw a GUI::Exception; the .what() value of the
                 * exception will be displayed to the user as an error message. */
                std::function<void(GUI::Parameters params)> &setup,
                /** Called to resume the simulation for count periods. */
                std::function<void(unsigned int rounds)> &resume,
                /** Called when the user hits the stop button in the GUI. This should cancel the
                 * current simulation run loop. */
                std::function<void()> &stop
        );

        /// Destructor: tells GUI to quit and rejoins GUI thread
        ~GUIShim();

        /** Checks whether the GUI has generated any events and, if so, processes them.  If there
         * are no pending events, this returns immediately.
         */
        void checkEvents();

        /** If the GUI has queued events, processes them.  If there are no queued events, this waits
         * for an event to occur then processes it.
         */
        void waitEvents();

        /** Starts the GUI. */
        void start(int argc, char *argv[]);

        /** Resynchronizes the simulation's agents into the GUI and tells the GUI to update its
         * graph.
         */
        void sync();

    protected:
        /** Processes a single event from the event queue and call the appropriate handler. */
//        processEvent()

        std::shared_ptr<Simulation> sim_;
        GUI gui_;

};

}
