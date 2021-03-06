#include "creativity/Creativity.hpp"
#include "creativity/state/State.hpp"
#include "creativity/state/Storage.hpp"
#include "creativity/state/FileStorage.hpp"
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
        creativity.read<FileStorage>(args.input, FileStorage::Mode::READONLY, args.memory_xz, args.tmpdir);
    }
    catch (std::exception &e) {
        std::cerr << "Unable to read `" << args.input << "': " << e.what() << "\n\n";
        exit(1);
    }

    auto stpair = creativity.storage();
    auto &st = *stpair.first;
    auto disp_size = st.size();
    if (disp_size > 0) disp_size--;

    std::cout << "Initial settings:\n=================\n" << std::left << std::setprecision(args.output_precision);
#define PRINT_FIELD(N, V) do { std::cout << std::setw(30) << N << V << "\n"; } while (0)
#define PRINT_SETTING(S) PRINT_FIELD(#S, creativity.parameters.S)
#define PRINT_SETTING_ARRAY(S) do { std::cout << std::setw(30) << #S << "["; \
    for (size_t i = 0; i < creativity.parameters.S.size(); i++) { \
        if (i > 0) std::cout << ", "; \
        std::cout << creativity.parameters.S[i]; \
    } \
    std::cout << "]\n"; } while (0)
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

    const auto &policies = creativity.parameters.policy;
    std::cout << std::setw(30) << "policy";
    bool comma = false;
    if (policies.publicSharing()) { if (comma) std::cout << ", "; comma = true; std::cout << "public sharing"; }
    if (policies.publicVoting())  { if (comma) std::cout << ", "; comma = true; std::cout << "public voting"; }
    if (policies.catchPirates())  { if (comma) std::cout << ", "; comma = true; std::cout << "catch pirates"; }
    if (policies.unknown())       { if (comma) std::cout << ", "; comma = true; std::cout << "(unknown policy)"; }
    if (!policies) std::cout << "no policy response";
    std::cout << "\n";
    PRINT_SETTING(policy_begins);

    PRINT_SETTING(policy_public_sharing_tax);

    PRINT_SETTING(policy_catch_tax);
    PRINT_SETTING_ARRAY(policy_catch_fine);
    PRINT_SETTING_ARRAY(policy_catch_mu);
    PRINT_SETTING_ARRAY(policy_catch_sigma);

    PRINT_SETTING(prior_scale);
    PRINT_SETTING(prior_scale_burnin);
    PRINT_SETTING(prior_scale_piracy);
    PRINT_SETTING(prior_scale_policy);
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
        std::cout << "\n\nArguments to re-run: ./creativity-cli --periods " << disp_size;
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

        std::list<std::string> policies;
        if (creativity.parameters.policy.publicSharing()) policies.push_back("public-sharing");
        if (creativity.parameters.policy.publicVoting())  policies.push_back("public-voting");
        if (creativity.parameters.policy.catchPirates())  policies.push_back("catch-pirates");
        if (creativity.parameters.policy.unknown())       policies.push_back("unknown-policy");
        if (!policies.empty()) {
            std::cout << " --policy ";
            bool first = true;
            for (const auto &p : policies) {
                if (first) first = false;
                else std::cout << ",";
                std::cout << p;
            }

            ADD_SAME(policy_begins);
        }

        if (creativity.parameters.policy.publicSharing())
            ADD("--public-sharing-tax", policy_public_sharing_tax);

        if (creativity.parameters.policy.publicVoting()) {
            ADD("--public-voting-tax", policy_public_voting_tax);
            ADD("--public-voting-votes", policy_public_voting_votes);
        }

        if (creativity.parameters.policy.catchPirates()) {
            ADD("--catch-tax", policy_catch_tax);
            ADD("--catch-fine-any", policy_catch_fine[0]);
            ADD("--catch-fine-const", policy_catch_fine[1]);
            ADD("--catch-fine-lin", policy_catch_fine[2]);
            ADD("--catch-fine-sq", policy_catch_fine[3]);
            ADD("--catch-mu-const", policy_catch_mu[0]);
            ADD("--catch-mu-lin", policy_catch_mu[1]);
            ADD("--catch-mu-sq", policy_catch_mu[2]);
            ADD("--catch-sigma-const", policy_catch_sigma[0]);
            ADD("--catch-sigma-lin", policy_catch_sigma[1]);
            ADD("--catch-sigma-sq", policy_catch_sigma[2]);
        }
        ADD_SAME(prior_scale);
        ADD_SAME(prior_scale_burnin);
        ADD_SAME(prior_scale_piracy);
        if (!policies.empty()) ADD_SAME(prior_scale_policy);
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


    std::cout << "\n\nSimulation details:\n===================\n";
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
             did_policy = creativity.parameters.policy_begins == 0;
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
            if (not did_policy and s->t >= creativity.parameters.policy_begins) {
                if (s->t == creativity.parameters.policy_begins) std::cout <<
                    "  =======================================  POLICY BEGINS  ========================================\n";
                else std::cout <<
                    "  ===================================  POLICY BEGINS " << std::setw(8) << std::left << std::string("(t=" + std::to_string(
                        creativity.parameters.policy_begins) + ")") <<           "  ===================================\n";
                did_policy = true;
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

