#include "creativity/Creativity.hpp"
#include "creativity/state/State.hpp"
#include "creativity/state/Storage.hpp"
#include "creativity/cmdargs/Info.hpp"
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <memory>
#include <string>
#include <map>
#include <utility>
#include <iostream>
#include <iomanip>

using namespace creativity;
using namespace creativity::state;
using namespace eris;
using namespace Eigen;

int main(int argc, char *argv[]) {
    cmdargs::Info args;
    args.parse(argc, argv);

    auto creativity = Creativity::create();
    // Filename input
    try {
        creativity->fileRead(args.input);
    }
    catch (std::exception &e) {
        std::cerr << "Unable to read `" << args.input << "': " << e.what() << "\n\n";
        exit(1);
    }

    std::cout << "Initial settings:\n=================\n" << std::left << std::setprecision(args.output_precision);
#define PRINT_FIELD(N, V) do { std::cout << std::setw(30) << N << V << "\n"; } while (0)
#define PRINT_SETTING(S) PRINT_FIELD(#S, creativity->parameters.S)
    PRINT_SETTING(readers);
    PRINT_SETTING(dimensions);
    PRINT_FIELD("densityFromBoundary()", creativity->densityFromBoundary());
    PRINT_SETTING(boundary);
    PRINT_SETTING(book_distance_sd);
    PRINT_SETTING(book_quality_sd);
    PRINT_SETTING(reader_step_sd);
    PRINT_SETTING(reader_creation_shape);
    PRINT_SETTING(reader_creation_scale_min);
    PRINT_SETTING(reader_creation_scale_max);
    PRINT_SETTING(creation_time);
    PRINT_SETTING(cost_fixed);
    PRINT_SETTING(cost_unit);
    PRINT_SETTING(cost_piracy);
    PRINT_SETTING(income);
    PRINT_SETTING(piracy_begins);
    PRINT_SETTING(piracy_link_proportion);
    PRINT_SETTING(prior_scale);
    PRINT_SETTING(prior_scale_piracy);
    PRINT_SETTING(prior_scale_burnin);
    PRINT_SETTING(burnin_periods);
    PRINT_SETTING(initial.prob_write);
    PRINT_SETTING(initial.q_min);
    PRINT_SETTING(initial.q_max);
    PRINT_SETTING(initial.p_min);
    PRINT_SETTING(initial.p_max);
    PRINT_SETTING(initial.prob_keep);
    PRINT_SETTING(initial.keep_price);
    PRINT_SETTING(initial.belief_threshold);
#undef PRINT_SETTING

    auto stpair = creativity->storage();
    auto &st = *stpair.first; 
    std::cout << "\n\nSimulation details:\n===================\n";
    auto disp_size = st.size();
    if (disp_size > 0) disp_size--;
    std::cout << "Periods: " << disp_size << "\n";
    if (st.size() > 0) {
        std::cout << "Period summary:\n";
        auto h1 =  "             Total                                         Average                                \n";
        auto h2 =  "       -----------------   -----------------------------------------------------------------------\n";
        auto h3 =  "   t   Rdrs  Books  BNew   Net.Util.  NewBuys  NewPiracy  NewPublic  TotBuys  TotPiracy  TotPublic\n";
        auto h4 =  "  ---  ----  -----  ----   ---------  -------  ---------  ---------  -------  ---------  ---------\n";
        auto fmt = "  %3d  %4d  %5d  %4d   %+9.3f  %7.3f  %9.3f  %9.3f  %7.2f  %9.2f  %9.2f\n";

        std::cout << h1 << h2;

        for (auto &s : st) {
            if (s->t == 0) continue;
            if (s->t % 25 == 1) { if (s->t > 1) std::cout << h4; std::cout << h3 << h4; }
            if (s->t == creativity->parameters.piracy_begins)
                std::cout <<
                   "  =======================================  PIRACY BEGINS  ========================================\n";
            if (s->t == creativity->parameters.public_sharing_begins)
                std::cout <<
                   "  ===================================  PUBLIC SHARING BEGINS  ====================================\n";

            double net_u = 0;
            unsigned long books_new = 0, bought = 0, publicb = 0, pirated = 0, bought_new = 0, pirated_new = 0, public_new = 0;
            for (auto &rp : s->readers) {
                auto &r = rp.second;
                net_u += r.u - 1000;

                bought += r.library_purchased;
                bought_new += r.library_purchased_new;
                publicb += r.library_public;
                public_new += r.library_public_new;
                pirated += r.library_pirated;
                pirated_new += r.library_pirated_new;
            }
            for (auto &bp : s->books) {
                auto &b = bp.second;
                if (s->t == b.created) books_new++;
            }

            double R = s->readers.size();
            printf(fmt,
                    s->t,
                    s->readers.size(),
                    s->books.size(),
                    books_new,
                    net_u / R,
                    bought_new / R,
                    pirated_new / R,
                    public_new / R,
                    bought / R,
                    pirated / R,
                    publicb / R);
        }
    }
    std::cout << "\n\n\n";
}

