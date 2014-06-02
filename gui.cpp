#include "creativity/GUI.hpp"
#include "creativity/Reader.hpp"
#include "creativity/Book.hpp"
#include "creativity/common.hpp"
#include "creativity/belief/Demand.hpp"
#include "creativity/belief/Profit.hpp"
#include "creativity/belief/ProfitStream.hpp"
#include "creativity/belief/Quality.hpp"
#include <eris/Eris.hpp>
#include <eris/Simulation.hpp>
#include <eris/Random.hpp>
#include <eris/intraopt/FixedIncome.hpp>
#include <functional>
#include <iostream>
#include <iomanip>
#include <Eigen/Core>

using namespace creativity;
using namespace eris;
using namespace Eigen;

int main(int argc, char *argv[1]) {
    ERIS_DBG("WOO");
    Eris<Simulation> sim;
    ERIS_DBG("WOO");
    MONEY = sim->create<Good::Continuous>();
    ERIS_DBG("WOO");
    sim->create<NEW_BOOKS_Cleaner>();

    ERIS_DBG("WOO");
    std::cerr << std::setprecision(16);
    std::cout << std::setprecision(16);

    ERIS_DBG("WOO");
    bool setup = false, stopped = true, step = false, quit = false;
    unsigned long num_readers = 1000;
    double book_sd = 0.5, quality_draw_sd = 1.0;
    double cost_fixed = 20, cost_unit = 1, income = 1000;
    unsigned long run_start = 0, run_end = 0;
    double speed_limit = 0;
    std::chrono::milliseconds redraw{50};

    ERIS_DBG("WOO");
#define NONNEG_DOUBLE(VAR, DESC) \
    case GUI::ParamType::VAR: \
        if (p.dbl < 0) \
            throw std::domain_error{DESC " `" + std::to_string(p.dbl) + "' is invalid"}; \
        VAR = p.dbl; \
        break;

    // Set up handlers for user actions in the GUI
    auto on_setup = [&](GUI::Parameter p) { // Setup
        ERIS_DBG("");
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
            NONNEG_DOUBLE(book_sd, "Book standard deviation");
            NONNEG_DOUBLE(quality_draw_sd, "Quality standard deviation");
            NONNEG_DOUBLE(cost_fixed, "Fixed cost");
            NONNEG_DOUBLE(cost_unit, "Unit cost");
            case GUI::ParamType::redraw:
                redraw = p.dur_ms;
                break;
            case GUI::ParamType::speed_limit:
                speed_limit = p.dbl;
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
        ERIS_DBG("done setup");
    };
    ERIS_DBG("WOO");
    auto on_run = [&](unsigned long periods) { // Run
        ERIS_DBG("run");
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
    ERIS_DBG("WOO");
    GUI gui(sim, on_setup, on_run, on_stop, on_resume, on_step, on_quit);

    try {
        gui.start(argc, argv);
    }
    catch (Glib::Error &e) {
        std::cerr << "Unable to start gui: " << e.what() << "\n";
        throw;
    }
    ERIS_DBG("WOO");

    auto &rng = eris::Random::rng();
    std::uniform_real_distribution<double> unif_01{0, 1};
    std::uniform_real_distribution<double> unif_pmb{-BOUNDARY, BOUNDARY};

    while (!setup) {
        gui.waitEvents();
        if (quit) return 0;
    }

    ERIS_DBG("Setting up readers");
    VectorXd demand_beta{8}; demand_beta << 0, -2, 0.5, -0.1, 0, 0, 0, 0;
    MatrixXd demand_V = MatrixXd::Identity(8, 8);
    double demand_s2 = 10, demand_n = 1;
    VectorXd profit_beta{5}; profit_beta << 0, 1, 0, 0, 0;
    MatrixXd profit_V = MatrixXd::Identity(5, 5);
    double profit_s2 = 10, profit_n = 1;
    VectorXd quality_beta{7}; quality_beta << 5, -1, 1, 0, 0, 0, 0.1;
    MatrixXd quality_V = MatrixXd::Identity(7, 7);
    double quality_s2 = 10, quality_n = 1;
    VectorXd stream_beta = VectorXd::Zero(10);
    MatrixXd stream_V = 1e10 * MatrixXd::Identity(10, 10);
    double stream_s2 = 1e20, stream_n = 0.0001;
    for (auto i = 0UL; i < num_readers; i++) {
        auto r = sim->create<Reader>(Position{unif_pmb(rng), unif_pmb(rng)},
                Position{-BOUNDARY,-BOUNDARY}, Position{BOUNDARY, BOUNDARY},
                belief::Demand{2, demand_beta, demand_s2, demand_V, demand_n},
                belief::Profit{2, profit_beta, profit_s2, profit_V, profit_n},
                belief::ProfitStream{stream_beta, stream_s2, stream_V, stream_n},
                belief::Quality{quality_beta, quality_s2, quality_V, quality_n},
                cost_fixed, cost_unit, income
                );
        r->writer_book_sd = book_sd;
    }
    ERIS_DBG("Done with readers");
    auto readers = sim->agents<Reader>();

    // Tell the GUI we're done with initialization
    gui.initialized();
    ERIS_DBG("here we go...!");

    gui.progress(run_end, 0);
    auto last_progress = std::chrono::high_resolution_clock::now();
    auto last_progress_t = sim->t();
    gui.redraw(true);

    // If redraw is 0, we use 
    constexpr auto zero_ms = std::chrono::milliseconds::zero();
    const auto never = std::chrono::time_point<std::chrono::high_resolution_clock>::max();
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
    }
    std::cerr << "running off the bottom of main\n";
}
