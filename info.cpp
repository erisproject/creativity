#include "creativity/Creativity.hpp"
#include "creativity/state/Storage.hpp"
#include "creativity/state/PsqlStorage.hpp"
#include <cstdio>

using namespace creativity;
using namespace creativity::state;
using namespace eris;
using namespace Eigen;

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " DATASOURCE -- print simulation summary information\n\n";
        std::cerr << "DATASOURCE can be a filename (typically a .crstate file), or a database URL,\n";
        std::cerr << "for example: 'postgresql://user:secret@localhost:5432/dbname?sslmode=require&creativity=123'\n\n\n";
        exit(1);
    }

    std::string source(argv[1]);
    auto creativity = Creativity::create();
    if (source.substr(0, 13) == "postgresql://" or source.substr(0, 11) == "postgres://") {
        try {
            creativity->pgsql(source, true /*read-only*/);

            PsqlStorage &pgsql = dynamic_cast<PsqlStorage&>(creativity->storage().first->backend());
            auto conn_locked = pgsql.connection();
            auto &conn = conn_locked.first;
            std::cout << "Connected to postgresql://";
            const char *username = conn.username();
            if (username) std::cout << username << "@";
            const char *hostname = conn.hostname();
            if (hostname) std::cout << hostname;
            const char *port = conn.port();
            if (port) std::cout << ":" << port;
            std::cout << "/" << conn.dbname() << "?creativity=" << pgsql.id << "\n";
        }
        catch (std::exception &e) {
            std::cerr << "Unable to connect to database `" << source << "': " << e.what() << "\n\n";
            exit(1);
        }
    }
    else {
        // Filename input
        try {
            creativity->fileRead(argv[1]);
        }
        catch (std::exception &e) {
            std::cerr << "Unable to read `" << argv[1] << "': " << e.what() << "\n\n";
            exit(1);
        }
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
        auto h1 =  "             Total                              Average                     \n";
        auto h2 =  "       -----------------   -------------------------------------------------\n";
        auto h3 =  "   t   Rdrs  Books  BNew   Net.Util.  NewBuys  NewPiracy  TotBuys  TotPiracy\n";
        auto h4 =  "  ---  ----  -----  ----   ---------  -------  ---------  -------  ---------\n";
        auto fmt = "  %3d  %4d  %5d  %4d   %+9.3f  %7.3f  %9.3f  %7.2f  %9.2f\n";

        std::cout << h1 << h2;

        for (auto &s : st) {
            if (s->t == 0) continue;
            if (s->t % 25 == 1) { if (s->t > 1) std::cout << h4; std::cout << h3 << h4; }
            if (s->t == creativity->parameters.piracy_begins)
                std::cout <<
                   "  ----------------------------  PIRACY BEGINS  -----------------------------\n";
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

