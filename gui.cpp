#include "creativity/GUI.hpp"
#include "creativity/Reader.hpp"
#include "creativity/Book.hpp"
#include "creativity/common.hpp"
#include <eris/Eris.hpp>
#include <eris/Simulation.hpp>
#include <eris/Random.hpp>
#include <eris/intraopt/FixedIncome.hpp>
#include <functional>
#include <iostream>

using namespace creativity;
using namespace eris;

int main(int argc, char *argv[1]) {
    Eris<Simulation> sim;
    MONEY = sim->create<Good::Continuous>();
    sim->create<NEW_BOOKS_Cleaner>();

    bool setup = false, stopped = false, step = false, quit = false;
    unsigned long num_readers = 1000;
    double prob_writer = 0.001;
    double book_sd = 0.5;
    unsigned long run_start = 0, run_end = 999;
    std::chrono::milliseconds speed_limit{50};
    std::chrono::milliseconds redraw{50};

    // Set up handlers for user actions in the GUI
    auto on_setup = [&](GUI::Parameter p) { // Setup
        switch (p.param) {
            case GUI::ParamType::begin:
                setup = false;
                break;
            case GUI::ParamType::dimensions:
                if (p.ul != 2)
                    throw std::domain_error{"Cannot yet handle dimensions ≠ 2"};
                // FIXME
                break;
            case GUI::ParamType::readers:
                if (p.ul < 1)
                    throw std::domain_error{"Must have at least one reader"};
                num_readers = p.ul;
                break;
            case GUI::ParamType::prob_writer:
                if (p.dbl < 0 || p.dbl > 1)
                    throw std::domain_error{"Book probability `" + std::to_string(p.dbl) + "' is invalid"};
                prob_writer = p.dbl;
                break;
            case GUI::ParamType::book_sd:
                if (p.dbl < 0)
                    throw std::domain_error{"Book standard deviation `" + std::to_string(p.dbl) + "' is invalid"};
                book_sd = p.dbl;
                break;
            case GUI::ParamType::redraw:
                redraw = p.dur_ms;
                break;
            case GUI::ParamType::speed_limit:
                speed_limit = p.dur_ms;
                break;
            case GUI::ParamType::threads:
                sim->maxThreads(p.ul);
                break;
            case GUI::ParamType::erred:
                break;
            case GUI::ParamType::finished:
                setup = true;
                break;
        }
    };
    auto on_run = [&](unsigned long periods) { // Run
        if (not setup)
            throw std::logic_error{"Event error: RUN before successful SETUP"};
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
    std::uniform_real_distribution<double> unif_pmb{-BOUNDARY, BOUNDARY};

    while (!setup) {
        gui.waitEvents();
        if (quit) return 0;
    }

    for (auto i = 0UL; i < num_readers; i++) {
        auto r = sim->create<Reader>(Position{unif_pmb(rng), unif_pmb(rng)},
                Position{-BOUNDARY,-BOUNDARY}, Position{BOUNDARY, BOUNDARY});
        r->writer_prob = prob_writer;
        r->writer_book_sd = book_sd;
        sim->create<intraopt::FixedIncome>(r, Bundle{{ MONEY, 1000 }});
    }
    auto readers = sim->agents<Reader>();

    gui.progress(run_end, 0);
    auto last_progress = std::chrono::high_resolution_clock::now();
    auto last_progress_t = sim->t();
    gui.redraw(true);

    // If redraw is 0, we use 
    constexpr auto zero_ms = std::chrono::milliseconds::zero();
    constexpr auto never = std::chrono::time_point<std::chrono::high_resolution_clock>::max();
    constexpr auto progress_freq = std::chrono::milliseconds{50};

    std::chrono::time_point<std::chrono::high_resolution_clock> start, end,
        next_sync{redraw == zero_ms ? never : std::chrono::high_resolution_clock::now() + redraw};

    std::chrono::time_point<std::chrono::high_resolution_clock> next_progress{
        std::chrono::high_resolution_clock::now() + progress_freq};

    while (not quit) {
        std::cerr << "asfafasdf\n";
        if (sim->t() < run_end) {
            // Tell the GUI we've started running.
            gui.running();
        }
        while (not quit and (step or (not stopped and sim->t() < run_end))) {
            start = std::chrono::high_resolution_clock::now();

            std::cerr << "running\n";
            sim->run();
            std::cerr << "done running\n";

            if (step) step = false;

            bool finished = stopped or sim->t() >= run_end;
            ERIS_DBG("");
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

            if (quit) break;

            // Make sure stopped hasn't changed:
            if (not finished and stopped) finished = true;

            // Only trigger a redraw and event check at most once every 50ms, or if we're done
            if (finished or end >= next_sync) {
                gui.redraw(true);
                next_sync = redraw == zero_ms ? never : end + redraw;
            }

            ERIS_DBG("");
            auto sleep = speed_limit - (end - start);
            if (not finished and sleep > zero_ms) {
                std::this_thread::sleep_for(sleep);
            }
            ERIS_DBG("");
        }
            ERIS_DBG("");

        // If the GUI told us to quit, just quit.
        if (quit) break;

        // Tell the GUI we finished
        gui.stopped(sim->t() >= run_end);

        std::cerr << "waiting for more\n";
        // Wait for the GUI to tell us to do something else
        gui.waitEvents();
        std::cerr << "got more\n";
    }
    std::cerr << "running off the bottom of main\n";
}
