#include "creativity/Creativity.hpp"
#include "creativity/CmdArgs.hpp"
#include "creativity/state/Storage.hpp"
#include "creativity/gui/GUI.hpp"
#include <eris/Simulation.hpp>
#include <eris/Random.hpp>
#include <Eigen/Core>
#include <functional>
#include <iostream>
#include <iomanip>
#include <cstddef>
#include <algorithm>
#include <chrono>
#include <ratio>
#include <stdexcept>
#include <string>
#include <utility>
#include <glibmm/error.h>
#include <glibmm/ustring.h>


using namespace creativity;
using namespace creativity::gui;
using namespace creativity::state;
using namespace eris;
using namespace Eigen;

int main(int argc, char *argv[1]) {
    Eigen::initParallel();
    auto creativity = Creativity::create();

    // Handle any command line arguments.  These set things in creativity->set(), which form the
    // defaults for the GUI, and have a few extra leftovers that the GUI code handles itself.
    CmdArgs::GUI cmd(creativity->set());
    cmd.parse(argc, argv);

    std::cerr << std::setprecision(16);
    std::cout << std::setprecision(16);

    bool setup = false, stopped = true, step = false, quit = false;
    eris_time_t run_start = 0, run_end = 0;
    std::chrono::milliseconds sync_speed{50};
    bool save_to_file = false, load_from_file = false;
    unsigned int max_threads = cmd.threads;

    // Set up handlers for user actions in the GUI
    auto on_configure = [&](GUI::Parameter p) { // Parameter
        switch (p.param) {
            case GUI::ParamType::seed:
                if (setup or load_from_file) throw std::runtime_error("Cannot change seed after initial setup");
                eris::Random::seed(p.ul);
                break;
            case GUI::ParamType::load:
                if (setup) throw std::runtime_error("Cannot load after initial setup");
                if (save_to_file) throw std::runtime_error("Error: cannot load and save at the same time");
                creativity->fileRead(*reinterpret_cast<std::string*>(p.ptr));
                if (creativity->storage().first->size() == 0)
                    throw std::runtime_error("Unable to load file: file has no states");
                if (creativity->parameters.dimensions < 1)
                    throw std::runtime_error("Unable to load file: invalid dimensions < 1");
                if (creativity->parameters.boundary <= 0)
                    throw std::runtime_error("Unable to load file: file has invalid non-positive boundary value");
                load_from_file = true;
                break;
            case GUI::ParamType::save_as:
                // TODO: this could actually be handled when already setup: it should be possible
                // to copy the current Storage object into the new FileStorage.
                if (setup) throw std::runtime_error("Cannot change file after initial setup");
                if (load_from_file) throw std::runtime_error("Error: cannot load and save at the same time");
                creativity->fileWrite(*reinterpret_cast<std::string*>(p.ptr));
                save_to_file = true;
                break;
            case GUI::ParamType::threads:
                // This is the only setting that *can* be changed after the initial setup.  This
                // will throw if currently running, but that's okay: the GUI isn't allowed to send
                // it if the simulation is currently running.
                if (creativity->sim) creativity->sim->maxThreads(p.ul);
                else max_threads = p.ul;
                break;
        }
    };

    auto on_initialize = [&]() {
        if (load_from_file)
            throw std::logic_error("Cannot initialize a new simulation after loading one from a file");
        if (setup)
            throw std::logic_error("Cannot initialize: initialization already done!");
        if (creativity->parameters.dimensions < 1) throw std::domain_error(u8"Invalid simulation file: dimensions < 1");
        // setup() will check the other parameters for validity
        creativity->setup();
        creativity->sim->maxThreads(max_threads);
        setup = true;
    };

    auto on_run = [&](eris_time_t periods) { // Run
        if (not setup)
            throw std::logic_error("Event error: RUN before successful SETUP");
        if (load_from_file)
            throw std::logic_error("Error: simulation state is readonly");
        run_start = creativity->sim->t();
        run_end = run_start + periods;
        stopped = false;
    };
    auto on_stop = [&]() { stopped = true; };
    auto on_step = [&]() { step = true; };
    auto on_resume = [&]() { stopped = false; };
    auto on_quit = [&]() { quit = true; };

    GUI gui(creativity, on_configure, on_initialize, on_run, on_stop, on_resume, on_step, on_quit);

    try {
        gui.start(cmd);
    }
    catch (Glib::Error &e) {
        std::cerr << "Unable to start gui: " << e.what() << "\n";
        throw;
    }

    while (!setup) {
        gui.waitEvents();
        if (quit) return 0;
    }
    if (load_from_file) {
        // We're in read-only mode, which means we do nothing but wait for the GUI to quit.
        gui.initialized();
        size_t num_states = creativity->storage().first->size() - 1; // -1 because we don't count the initial, t=0 setup state
        gui.progress(num_states, num_states, 0);
        gui.newStates(0);
        while (not quit) {
            gui.waitEvents();
        }
        return 0;
    }

    auto sim = creativity->sim;

    // Copy the initial state into the storage object
    creativity->storage().first->emplace_back(sim);

    // Tell the GUI we're done with initialization so that it can disable the various setup elements
    gui.initialized();

    gui.progress(0, run_end, 0);
    auto last_progress = std::chrono::high_resolution_clock::now();
    auto last_progress_t = sim->t();

    constexpr auto progress_freq = std::chrono::milliseconds{50};

    std::chrono::time_point<std::chrono::high_resolution_clock> start, end,
        next_sync{std::chrono::high_resolution_clock::now() + sync_speed};

    std::chrono::time_point<std::chrono::high_resolution_clock> next_progress{
        std::chrono::high_resolution_clock::now() + progress_freq};

    while (not quit) {
        if (sim->t() < run_end) {
            // Tell the GUI we've started running.
            gui.running();
        }
        while (not quit and (step or (not stopped and sim->t() < run_end))) {
            start = std::chrono::high_resolution_clock::now();

            creativity->run();

            if (step) step = false;

            bool finished = stopped or sim->t() >= run_end;
            end = std::chrono::high_resolution_clock::now();

            if (finished or end >= next_progress) {
                auto now_t = sim->t();
                double speed = (double) (now_t - last_progress_t) / std::chrono::duration<double>{end - last_progress}.count();
                gui.progress(now_t, std::max(now_t, run_end), speed);
                last_progress = end;
                last_progress_t = now_t;
                gui.checkEvents();
                next_progress = end + progress_freq;
            }

            if (quit) break;

            // Make sure stopped hasn't changed:
            if (not finished and stopped) finished = true;

            // Only tell the GUI about new states at most once every 50ms, or if we're done
            if (finished or end >= next_sync) {
                gui.newStates();
                next_sync = end + sync_speed;
            }
        }

        // If the GUI told us to quit, just quit.
        if (quit) break;

        // Tell the GUI we finished
        gui.stopped(sim->t() < run_end);

        // Wait for the GUI to tell us to do something else
        gui.waitEvents();

        // Update this so that the speed is right (otherwise however long waitEvents() waits for the
        // user waits would be included in the calculation).
        last_progress = std::chrono::high_resolution_clock::now();
    }

    creativity->storage().first->flush();
}
