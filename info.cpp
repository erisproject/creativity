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
#include <regex>

using namespace creativity;
using namespace creativity::state;
using namespace eris;
using namespace Eigen;

int main(int argc, char *argv[]) {
    cmdargs::Info args;
    args.parse(argc, argv);

    Creativity creativity;
    // Filename input
    try {
        creativity.fileRead(args.input, args.memory_xz or args.memory, args.memory);
    }
    catch (std::exception &e) {
        std::cerr << "Unable to read `" << args.input << "': " << e.what() << "\n\n";
        exit(1);
    }

    std::cout << "Initial settings:\n=================\n" << std::left << std::setprecision(args.output_precision);
#define PRINT_FIELD(N, V) do { std::cout << std::setw(30) << N << V << "\n"; } while (0)
#define PRINT_SETTING(S) PRINT_FIELD(#S, creativity.parameters.S)
    PRINT_SETTING(readers);
    PRINT_SETTING(dimensions);
    PRINT_FIELD("densityFromBoundary()", creativity.densityFromBoundary());
    PRINT_SETTING(boundary);
    PRINT_SETTING(book_distance_mean);
    PRINT_SETTING(book_quality_sd);
    PRINT_SETTING(reader_step_mean);
    PRINT_SETTING(reader_creation_shape);
    PRINT_SETTING(reader_creation_scale_min);
    PRINT_SETTING(reader_creation_scale_range);
    PRINT_SETTING(creation_time);
    PRINT_SETTING(creation_fixed);
    PRINT_SETTING(cost_market);
    PRINT_SETTING(cost_unit);
    PRINT_SETTING(cost_piracy);
    PRINT_SETTING(income);
    PRINT_SETTING(piracy_begins);
    PRINT_SETTING(piracy_link_proportion);
    PRINT_SETTING(public_sharing_begins);
    PRINT_SETTING(public_sharing_tax);
    PRINT_SETTING(prior_scale);
    PRINT_SETTING(prior_scale_burnin);
    PRINT_SETTING(prior_scale_piracy);
    PRINT_SETTING(prior_scale_public_sharing);
    PRINT_SETTING(burnin_periods);
    PRINT_SETTING(prediction_draws);
    PRINT_SETTING(initial.prob_write);
    PRINT_SETTING(initial.l_min);
    PRINT_SETTING(initial.l_range);
    PRINT_SETTING(initial.p_min);
    PRINT_SETTING(initial.p_range);
    PRINT_SETTING(initial.prob_keep);
    PRINT_SETTING(initial.keep_price);
    PRINT_SETTING(initial.belief_threshold);
#undef PRINT_SETTING

    if (args.show_cli_args) {
        CreativitySettings defaults;
        auto saveprec = std::cout.precision();
        std::cout.precision(std::numeric_limits<double>::max_digits10);
        std::cout << "\n\nArguments to re-run: ./creativity-cli";
#define ADD(ARG, PARAM) if (creativity.parameters.PARAM != defaults.PARAM) { \
    std::string arg(ARG); \
    std::replace(arg.begin(), arg.end(), '_', '-'); \
    std::replace(arg.begin(), arg.end(), '.', '-'); \
    std::cout << " " << arg << " " << creativity.parameters.PARAM; \
}
#define ADD_SAME(PARAM) ADD("--" #PARAM, PARAM)
        ADD_SAME(readers);
        ADD_SAME(dimensions);
        if (creativity.parameters.boundary != defaults.boundary) std::cout << " --density " <<
            creativity.densityFromBoundary(creativity.parameters.readers, creativity.parameters.dimensions, creativity.parameters.boundary);
        ADD_SAME(book_distance_mean);
        ADD_SAME(book_quality_sd);
        ADD_SAME(reader_step_mean);
        ADD_SAME(reader_creation_shape);
        ADD_SAME(reader_creation_scale_min);
        ADD_SAME(reader_creation_scale_range);
        ADD_SAME(creation_time);
        ADD_SAME(creation_fixed);
        ADD_SAME(cost_market);
        ADD_SAME(cost_unit);
        ADD_SAME(cost_piracy);
        ADD_SAME(income);
        ADD_SAME(piracy_begins);
        ADD_SAME(piracy_link_proportion);
        ADD_SAME(public_sharing_begins);
        ADD_SAME(public_sharing_tax);
        ADD_SAME(prior_scale);
        ADD_SAME(prior_scale_burnin);
        ADD_SAME(prior_scale_piracy);
        ADD_SAME(prior_scale_public_sharing);
        ADD_SAME(burnin_periods);
        ADD_SAME(prediction_draws);
        ADD_SAME(initial.prob_write);
        ADD("--initial-effort-min", initial.l_min);
        ADD("--initial-effort-range", initial.l_range);
        ADD("--initial-price-min", initial.p_min);
        ADD("--initial-price-range", initial.p_range);
        ADD_SAME(initial.prob_keep);
        ADD_SAME(initial.keep_price);
        ADD("--belief-threshold", initial.belief_threshold);

        std::smatch smatch;
        if (regex_search(args.input, smatch, std::regex("-(\\d+)\\.(?:crstate\\b|\\w+$)"))) {
            std::cout << " --seed " << smatch[1] << "\n";
        }
        else  {
            std::cout << "\n";
            std::cerr << "(could not determine --seed value from input filename `" << args.input << "')\n";
        }

        std::cout.precision(saveprec);
    }


    auto stpair = creativity.storage();
    auto &st = *stpair.first; 
    std::cout << "\n\nSimulation details:\n===================\n";
    auto disp_size = st.size();
    if (disp_size > 0) disp_size--;
    std::cout << "Periods: " << disp_size << "\n";
    if (st.size() > 0 and args.thin_periods > 0) {
        std::cout << "Period summary";
        if (args.thin_periods > 1) std::cout << " (every " << args.thin_periods <<
            (args.thin_periods >= 11 and args.thin_periods <= 13 ? "th" : args.thin_periods % 10 == 1 ? "st" :
             args.thin_periods % 10 == 2 ? "nd" : args.thin_periods % 10 == 3 ? "rd" : "th")
                << " period)";
        std::cout << ":\n";

        auto h1 =  "             Total                                         Average                                \n";
        auto h2 =  "       -----------------   -----------------------------------------------------------------------\n";
        auto h3 =  "   t   Rdrs  Books  BNew   Net.Util.  NewBuys  NewPiracy  NewPublic  TotBuys  TotPiracy  TotPublic\n";
        auto h4 =  "  ---  ----  -----  ----   ---------  -------  ---------  ---------  -------  ---------  ---------\n";
        auto fmt = "  %3d  %4d  %5d  %4d   %+9.3f  %7.3f  %9.3f  %9.3f  %7.2f  %9.2f  %9.2f\n";

        std::cout << h1 << h2 << h3 << h4;

        unsigned count = 0;
        // These are 0 for disabled, so pretend we already did them if 0:
        bool did_piracy = creativity.parameters.piracy_begins == 0,
             did_public_sharing = creativity.parameters.public_sharing_begins == 0;
        for (unsigned t = (args.thin_periods < st.size() ? args.thin_periods : st.size()-1); t < st.size(); t += args.thin_periods) {
            auto s = st[t];
            if (++count % 35 == 1 and count > 1) { std::cout << h4 << h3 << h4; }
            if (not did_piracy and s->t >= creativity.parameters.piracy_begins) {
                if (s->t == creativity.parameters.piracy_begins) std::cout <<
                   "  =======================================  PIRACY BEGINS  ========================================\n";
                else std::cout <<
                   "  ===================================  PIRACY BEGINS " << std::setw(8) << std::left << std::string("(t=" + std::to_string(
                    creativity.parameters.piracy_begins) + ")") <<            "  ===================================\n";
                did_piracy = true;
            }
            if (not did_public_sharing and s->t >= creativity.parameters.public_sharing_begins) {
                if (s->t == creativity.parameters.public_sharing_begins) std::cout <<
                   "  ===================================  PUBLIC SHARING BEGINS  ====================================\n";
                else std::cout <<
                   "  ===============================  PUBLIC SHARING BEGINS " << std::setw(8) << std::left << std::string("(t=" + std::to_string(
                    creativity.parameters.public_sharing_begins) + ")") <<         "  ===============================\n";
                did_public_sharing = true;
            }

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

            // Mess around with t if we're near but not at the end so that we always include the
            // very last period.
            if (t != st.size()-1 and t + args.thin_periods >= st.size()) {
                t = st.size()-1-args.thin_periods;
            }
        }
    }
    std::cout << "\n\n\n";
}

