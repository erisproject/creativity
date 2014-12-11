#include "creativity/Creativity.hpp"
#include "creativity/state/Storage.hpp"
#include <cstdio>

using namespace creativity;
using namespace creativity::state;
using namespace eris;
using namespace Eigen;

int main(int argc, char *argv[]) {
    if (argc <= 1) {
        std::cerr << "Usage: " << argv[0] << " FILENAME.crstate -- print summary information about a .crstate file\n\n";
        exit(1);
    }

    auto creativity = Creativity::create();
    try {
        creativity->fileRead(argv[1]);
    }
    catch (std::exception &e) {
        std::cerr << "Unable to read `" << argv[1] << "': " << e.what() << "\n\n";
        exit(1);
    }

    std::cout << "Initial settings:\n=================\n";
#define PRINT_SETTING(S) do { printf("%-30s", #S); std::cout << creativity->parameters.S << "\n"; } while (0)
    PRINT_SETTING(readers);
    PRINT_SETTING(dimensions);
    PRINT_SETTING(boundary);
    PRINT_SETTING(book_distance_sd);
    PRINT_SETTING(book_quality_sd);
    PRINT_SETTING(reader_step_sd);
    PRINT_SETTING(reader_creation_shape);
    PRINT_SETTING(reader_creation_scale_min);
    PRINT_SETTING(reader_creation_scale_max);
    PRINT_SETTING(cost_fixed);
    PRINT_SETTING(cost_unit);
    PRINT_SETTING(cost_piracy);
    PRINT_SETTING(income);
    PRINT_SETTING(piracy_begins);
    PRINT_SETTING(piracy_link_proportion);
    PRINT_SETTING(prior_weight);
    PRINT_SETTING(prior_weight_piracy);
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
    std::cout << "Periods: " << st.size()-1 << "\n";
    if (st.size() > 0) {
        std::cout << "Period summary:\n";
        printf(    "             Total                              Average                     \n");
        printf(    "       -----------------   -------------------------------------------------\n");
        printf(    "   t   Rdrs  Books  BNew   Net.Util.  NewBuys  NewPiracy  TotBuys  TotPiracy\n");
        printf(    "  ---  ----  -----  ----   ---------  -------  ---------  -------  ---------\n");
        auto fmt = "  %3d  %4d  %5d  %4d   %+9.3f  %7.3f  %9.3f  %7.2f  %9.2f\n";

        for (auto &s : st) {
            if (s->t == 0) continue;
            if (s->t == creativity->parameters.piracy_begins)
                std::cout <<
                   "  ------------------------------  PIRACY BEGINS  -------------------------------\n";
            double net_u = 0;
            unsigned long books_new = 0, bought = 0, pirated = 0, bought_new = 0, pirated_new = 0;
            for (auto &rp : s->readers) {
                auto &r = rp.second;
                net_u += r.u - 1000;
                bought += r.library_purchased.size();
                pirated += r.library_pirated.size();
                bought_new += r.new_purchased.size();
                pirated_new += r.new_pirated.size();
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
                    bought / R,
                    pirated / R);
        }
    }
    std::cout << "\n\n\n";
}

