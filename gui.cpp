#include "creativity/GUI.hpp"
#include "creativity/Reader.hpp"
#include "creativity/Book.hpp"
#include <eris/Eris.hpp>
#include <eris/Simulation.hpp>
#include <eris/Random.hpp>
#include <functional>
#include <iostream>

using namespace creativity;
using namespace std::placeholders;

int main(int argc, char *argv[1]) {
    Eris<Simulation> sim;

    bool setup = false, stopped = false, step = false, quit = false;
    unsigned long num_readers = 1000;
    double prob_writer = 0.001;
    double book_sd = 0.5;
    unsigned long run_start = 0, run_end = 999;
    std::chrono::milliseconds speed_limit{50};
    std::chrono::milliseconds redraw{50};

    // Set up handlers for user actions in the GUI
    auto on_setup = [&](GUI::Parameters p) { // Setup
        if (p.dimensions != 2)
            throw std::domain_error{"Cannot yet handle dimensions != 2"};
        if (p.readers < 1)
            throw std::domain_error{"Must have at least one reader"};
        if (p.prob_writer < 0 || p.prob_writer > 1)
            throw std::domain_error{"Book probability `" + std::to_string(p.prob_writer) + "' is invalid"};
        if (p.book_sd < 0)
            throw std::domain_error{"Book standard deviation `" + std::to_string(p.book_sd) + "' is invalid"};
        num_readers = p.readers;
        prob_writer = p.prob_writer;
        book_sd = p.book_sd;
        speed_limit = p.speed_limit;
        redraw = p.redraw;

        setup = true;
    };
    auto on_run = [&](unsigned long periods) { // Run
        if (not setup)
            throw std::logic_error{"Event error: RUN before SETUP"};
        run_start = sim->t();
        run_end = run_start + periods;
        stopped = false;
    };
    auto on_stop = [&]() { stopped = true; };
    auto on_step = [&]() { step = true; };
    auto on_resume = [&]() { stopped = false; };
    auto on_quit = [&]() { quit = true; };
    GUI gui(sim, on_setup, on_run, on_stop, on_resume, on_step, on_quit);

    try {
        gui.start(argc, argv);
    }
    catch (Glib::Error &e) {
        std::cerr << "Unable to start gui: " << e.what() << "\n";
        throw;
    }

    auto &rng = eris::Random::rng();
    std::uniform_real_distribution<double> unif_01{0, 1};
    std::uniform_real_distribution<double> unif_pm10{-10, 10};

    while (!setup) {
        gui.waitEvents();
        if (quit) return 0;
    }

    for (auto i = 0UL; i < num_readers; i++) {
        auto r = sim->create<Reader>(Position{unif_pm10(rng), unif_pm10(rng)});
        r->writer_prob = prob_writer;
        r->writer_book_sd = book_sd;
    }
    auto readers = sim->agents<Reader>();

    gui.progress(run_end, 0);
    auto last_progress = std::chrono::high_resolution_clock::now();
    auto last_progress_t = sim->t();
    gui.redraw();

    // If redraw is 0, we use 
    constexpr auto zero_ms = std::chrono::milliseconds::zero();
    constexpr auto never = std::chrono::time_point<std::chrono::high_resolution_clock>::max();
    constexpr auto progress_freq = std::chrono::milliseconds{50};

    std::chrono::time_point<std::chrono::high_resolution_clock> start, end,
        next_sync{redraw == zero_ms ? never : std::chrono::high_resolution_clock::now() + redraw};

    std::chrono::time_point<std::chrono::high_resolution_clock> next_progress{
        std::chrono::high_resolution_clock::now() + progress_freq};

    while (not quit) {
        if (sim->t() < run_end) {
            // Tell the GUI we've started running.
            gui.running();
        }
        while (not quit and (step or (not stopped and sim->t() < run_end))) {
            start = std::chrono::high_resolution_clock::now();

            sim->run();

            if (step) step = false;

            bool finished = stopped or sim->t() >= run_end;
            end = std::chrono::high_resolution_clock::now();
            // Only update the progress and check events every 50ms
            if (finished or end >= next_progress) {
                auto now_t = sim->t();
                double speed = (double) (now_t - last_progress_t) / std::chrono::duration<double>{end - last_progress}.count();
                gui.progress(run_end, speed);
                last_progress = end;
                last_progress_t = now_t;
                gui.checkEvents();
                next_progress = end + progress_freq;
            }

            // Only trigger a redraw and event check at most once every 50ms, or if we're done
            if (finished or end >= next_sync) {
                gui.redraw();
                next_sync = redraw == zero_ms ? never : end + redraw;
            }

            auto sleep = speed_limit - (end - start);
            if (not finished and sleep > zero_ms) {
                std::cerr << "sleep for " << sleep.count() << "\n";
                std::this_thread::sleep_for(sleep);
            }
        }

        // If the GUI told us to quit, just quit.
        if (quit) break;

        // Tell the GUI we finished
        gui.stopped(sim->t() >= run_end);

        // Wait for the GUI to tell us to do something else
        gui.waitEvents();
    }
}
