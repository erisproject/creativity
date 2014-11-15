#include "creativity/Creativity.hpp"
#include "creativity/gui/GUI.hpp"
#include <eris/Simulation.hpp>
#include <eris/Random.hpp>
#include <functional>
#include <iostream>
#include <iomanip>
#include <Eigen/Core>

using namespace creativity;
using namespace creativity::gui;
using namespace creativity::state;
using namespace eris;
using namespace Eigen;

int main(int argc, char *argv[1]) {
    Eigen::initParallel();
    auto creativity = Creativity::create();

    std::cerr << std::setprecision(16);
    std::cout << std::setprecision(16);

    bool setup = false, stopped = true, step = false, quit = false;
    unsigned long run_start = 0, run_end = 0;
    double speed_limit = 0;
    std::chrono::milliseconds sync_speed{50};
    bool save_to_file = false, load_from_file = false;
    unsigned int max_threads = 0;

#define NONNEG_DOUBLE(TYPE, VAR, DESC) \
    case GUI::ParamType::TYPE: \
        if (setup) throw std::runtime_error("Cannot change " DESC " after initial setup"); \
        if (p.dbl < 0) \
            throw std::domain_error{DESC " `" + std::to_string(p.dbl) + "' is invalid"}; \
        VAR = p.dbl; \
        break;

    // Set up handlers for user actions in the GUI
    auto on_setup = [&](GUI::Parameter p) { // Setup
        switch (p.param) {
            case GUI::ParamType::finished:
                creativity->setup();
                setup = true;
                break;
            case GUI::ParamType::begin:
            case GUI::ParamType::erred:
                break;
            case GUI::ParamType::dimensions:
                if (setup) throw std::runtime_error("Cannot change dimensions after initial setup");
                if (p.ul != 2)
                    throw std::domain_error{"Cannot yet handle dimensions ≠ 2"};
                // FIXME
                break;
            case GUI::ParamType::readers:
                if (setup) throw std::runtime_error("Cannot change readers after initial setup");
                if (p.ul < 1)
                    throw std::domain_error{"Must have at least one reader"};
                creativity->parameters.readers = p.ul;
                break;
            case GUI::ParamType::density:
                if (setup) throw std::runtime_error("Cannot change density after initial setup");
                if (p.dbl <= 0)
                    throw std::domain_error{"Density must be positive"};
                creativity->parameters.density = p.dbl;
                break;
            case GUI::ParamType::sharing_begins:
                if (setup) throw std::runtime_error("Cannot change sharing t after initial setup");
                creativity->parameters.sharing_begins = p.ul;
                break;
            case GUI::ParamType::seed:
                if (setup) throw std::runtime_error("Cannot change seed after initial setup");
                eris::Random::seed(p.ul);
                break;
            case GUI::ParamType::load:
                if (setup) throw std::runtime_error("Cannot load after initial setup");
                if (save_to_file) throw std::runtime_error("Error: setup specified both load and save file");
                creativity->fileRead(*reinterpret_cast<std::string*>(p.ptr));
                if (creativity->storage().first->size() == 0)
                    throw std::runtime_error("Unable to load file: file has no states");
                if (creativity->storage().first->dimensions() != 2)
                    throw std::runtime_error("Unable to load file: dimensions != 2");
                if (creativity->storage().first->boundary() <= 0)
                    throw std::runtime_error("Unable to load file: file has invalid non-positive boundary value");
                load_from_file = true;
                break;
            case GUI::ParamType::save_as:
                // FIXME: this could actually be handled when already setup: it should be possible
                // to copy the current Storage object into the new FileStorage.
                if (setup) throw std::runtime_error("Cannot change file after initial setup");
                if (load_from_file) throw std::runtime_error("Error: setup specified both load and save file");
                creativity->fileWrite(*reinterpret_cast<std::string*>(p.ptr));
                save_to_file = true;
                break;
            NONNEG_DOUBLE(book_sd, creativity->parameters.book_distance_sd, "Book standard deviation");
            NONNEG_DOUBLE(quality_draw_sd, creativity->parameters.book_quality_sd, "Quality standard deviation");
            NONNEG_DOUBLE(cost_fixed, creativity->parameters.cost_fixed, "Fixed cost");
            NONNEG_DOUBLE(cost_unit, creativity->parameters.cost_unit, "Unit cost");
            case GUI::ParamType::threads:
                // This is the only setting that *can* be changed after the initial setup.  This
                // will throw if currently running, but that's okay: the GUI isn't allowed to send
                // it if the simulation is currently running.
                max_threads = p.ul;
                break;
        }
    };
#undef NONNEG_DOUBLE

    auto on_run = [&](unsigned long periods) { // Run
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

    GUI gui(creativity, on_setup, on_run, on_stop, on_resume, on_step, on_quit);

    try {
        gui.start(argc, argv);
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
        gui.newStates(0);
        while (not quit) {
            gui.waitEvents();
        }
        return 0;
    }

    auto sim = creativity->sim;

    // Copy the initial state into the storage object
    creativity->storage().first->emplace_back(sim);

    // Tell the GUI we're done with initialization
    gui.initialized();
    ERIS_DBG("here we go...!");

    gui.progress(0, run_end, 0);
    auto last_progress = std::chrono::high_resolution_clock::now();
    auto last_progress_t = sim->t();

    constexpr auto zero_ms = std::chrono::milliseconds::zero();
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

            ERIS_DBG("running");
            sim->run();
            ERIS_DBG("done running");

            {
                ERIS_DBG("Adding simulation state to storage...");
                creativity->storage().first->emplace_back(sim);
                ERIS_DBG("done");
            }

            if (step) step = false;

            bool finished = stopped or sim->t() >= run_end;
            end = std::chrono::high_resolution_clock::now();

            if (finished or end >= next_progress) {
                auto now_t = sim->t();
                double speed = (double) (now_t - last_progress_t) / std::chrono::duration<double>{end - last_progress}.count();
                gui.progress(now_t, run_end, speed);
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

            if (not finished and speed_limit > 0) {

                auto sleep = std::chrono::duration<double>{1.0/speed_limit} - (end - start);
                if (sleep > zero_ms)
                    std::this_thread::sleep_for(sleep);
                // Else we're already slower than the speed limit
            }
        }

        // If the GUI told us to quit, just quit.
        if (quit) break;

        // Tell the GUI we finished
        gui.stopped(sim->t() < run_end);

        std::cerr << "waiting for more\n";
        // Wait for the GUI to tell us to do something else
        gui.waitEvents();
        std::cerr << "got more\n";

        // Update this so that the speed is right (otherwise however long waitEvents() waits for the
        // user waits would be included in the calculation).
        last_progress = std::chrono::high_resolution_clock::now();
    }

    creativity->storage().first->flush();

    std::cerr << "running off the bottom of main\n";
}
