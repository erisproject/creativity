#include "creativity/Creativity.hpp"
#include "creativity/Reader.hpp"
#include "creativity/Book.hpp"
#include "creativity/state/Storage.hpp"
#include "creativity/state/FileStorage.hpp"
#include "creativity/state/StorageBackend.hpp"
#include "creativity/cmdargs/CLI.hpp"
#include <eris/Simulation.hpp>
#include <Eigen/Core>
#include <cerrno>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <exception>
#include <queue>
#include <random>
#include <ratio>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#ifdef __linux__
extern "C" {
#include <sys/prctl.h>
}
#endif

using namespace creativity;
using namespace creativity::state;
using namespace eris;

namespace fs = boost::filesystem;

std::string random_filename(const std::string &basename) {
    std::random_device rd;
    boost::random::uniform_int_distribution<unsigned char> rand_index(0, 35);
    constexpr char map[37] = "0123456789abcdefghijklmnopqrstuvwxyz";

    std::ostringstream buf;
    buf << basename << ".partial.";
    for (int i = 0; i < 7; i++) {
        buf << map[rand_index(rd)];
    }
    return buf.str();
}

// Simple class used to format a duration on output
struct duration {
    unsigned days, hours, minutes, seconds, milliseconds;
    explicit duration(double s) :
        days{unsigned(s / 86400)},
        hours{unsigned((s - days*86400) / 3600)},
        minutes{unsigned((s - days*86400 - hours*3600) / 60)},
        seconds{unsigned(s - days*86400 - hours*3600 - minutes*60)},
        milliseconds{unsigned(std::lround(1000*(s - days*86400 - hours*3600 - minutes*60 - seconds)))}
    {}
    friend std::ostream& operator<<(std::ostream &out, const duration &d) {
        auto save = out.flags();
        if (d.days > 0)
            out << d.days << 'd';
        if (d.days > 0 or d.hours > 0)
            out << d.hours << 'h' << std::setfill('0') << std::setw(2);
        if (d.days > 0 or d.hours > 0 or d.minutes > 0)
            out << d.minutes << 'm' << std::setfill('0') << std::setw(2);
        out << d.seconds << '.' << std::setfill('0') << std::setw(3) << d.milliseconds << "s" << std::setfill(' ');
        out.flags(save);
        return out;
    }
};

int main(int argc, char *argv[]) {
    std::cerr << std::setprecision(16);
    std::cout << std::setprecision(16);

    auto cr_ptr = std::unique_ptr<Creativity>(new Creativity);
    auto &creativity = *cr_ptr;

    cmdargs::CLI args(creativity.set());
    args.parse(argc, argv);

    std::string args_out = args.output, tmp_out = args_out;
    bool overwrite = args.overwrite;

    // Filename output
    try {
        fs::path outpath(args_out);
        // If the user didn't specify --overwrite, make sure the file doesn't exist.  (It might
        // get created before we finish the simulation, but at least we can abort now in a
        // typical case).
        if (not overwrite and fs::exists(outpath)) {
            std::cerr << "Error: `" << args_out << "' already exists; specify a different file or add `--overwrite' option to overwrite\n";
            exit(1);
        }

        std::string basename = outpath.filename().string();
        std::string dirname = outpath.parent_path().string();
        if (dirname.empty()) dirname = ".";

        if (basename.empty()) throw std::runtime_error("Invalid output filename (" + args_out + ")");

        // Make sure the parent of the output file exists
        if (not fs::is_directory(dirname))
            throw std::runtime_error("Directory `" + dirname + "' does not exist or is not a directory");

        const auto &tmpdir = args.tmpdir;
        if (not tmpdir.empty() and not fs::is_directory(tmpdir))
            throw std::runtime_error("Error: --tmpdir `" + tmpdir + "' does not exist or is not a directory");

        creativity.write<FileStorage>(args_out, FileStorage::Mode::OVERWRITE, args.memory, tmpdir);
    }
    catch (std::exception &e) {
        std::cerr << "Unable to write to file: " << e.what() << "\n";
        exit(1);
    }

    creativity.setup();
    auto sim = creativity.sim;

    // Copy the initial state into the storage object
    creativity.storage().first->emplace_back(sim);

    if (args.threads > 0) Eigen::initParallel();
    sim->maxThreads(args.threads);

    constexpr size_t avg_times = 5;
    std::chrono::time_point<std::chrono::high_resolution_clock> started = std::chrono::high_resolution_clock::now();
    std::queue<std::chrono::time_point<std::chrono::high_resolution_clock>> times;
    times.push(started);
    unsigned int max_bnew_digits = 0;
    auto count_bnew = [](const Book &b) -> bool { return b.age() == 0; };

    unsigned int periods = args.periods;

    bool tty = isatty(fileno(stdout));
    if (args.quiet) std::cout << "Running simulation for " << periods << " periods.\n";
    while (sim->t() < periods) {
#ifdef __linux__
        // Update the process name to something like "crcli [43/300]" (the part before the [ gets
        // shorted if required--we aren't allowed to set a name longer than 15 characters).
        {
            std::string progress = "[" + std::to_string(sim->t()) + "/" + std::to_string(periods) + "]";
            std::string name = progress.size() <= 9 ? "crcli " + progress : "crcli" + progress;
            prctl(PR_SET_NAME, name.c_str());
        }
#endif
        creativity.run();

        if (not args.quiet) {
            // Calculate and show various status information, such as simulation speed, number of books,
            // etc.
            auto now = std::chrono::high_resolution_clock::now();
            std::ostringstream speed;
            speed << std::setprecision(6) << std::showpoint;
            speed << std::setw(7) << 1.0 / std::chrono::duration<double>(now - times.back()).count() << " Hz ";
            while (times.size() > avg_times) times.pop();
            if (times.size() > 1)
                speed << std::setw(7) << times.size() / std::chrono::duration<double>(now - times.front()).count() << " Hz[" << times.size() << "] ";

            auto bnew = sim->countGoods<Book>(count_bnew);
            max_bnew_digits = std::max(max_bnew_digits, bnew == 0 ? 1 : 1 + (unsigned) std::lround(std::floor(std::log10(bnew))));

            speed << std::setw(7) << sim->t() / std::chrono::duration<double>(now - started).count() << " Hz[A]";
            times.push(std::move(now));

            if (tty) std::cout << "\r";
            std::cout << "Running simulation [t=" << sim->t() << "; " <<
                (creativity.piracy() ? u8"Piracy ✔; " : u8"Piracy ✘; ") <<
                (creativity.policyActive() ? u8"Policy ✔; " : u8"Policy ✘; ") <<
                "R=" << sim->countAgents<Reader>() << "; B=" << sim->countGoods<Book>() << "; Bnew=" <<
                std::setw(max_bnew_digits) << bnew << "] " << speed.str();
            if (tty) std::cout << std::flush;
            else std::cout << std::endl;
        }
    }

    if (tty and not args.quiet) std::cout << std::endl;

    std::chrono::time_point<std::chrono::high_resolution_clock> finished_sim = std::chrono::high_resolution_clock::now();
    std::cout << "Simulation finished in ";
    std::cout << duration(std::chrono::duration<double>(finished_sim - started).count());
    std::cout << std::endl;

    bool need_endl = false;
    bool flush_waited = false;
    while (not creativity.storage().first->flush_for(250)) {
        if (args.quiet) {
            if (not flush_waited) {
                std::cout << "Waiting for output data to finish writing... " << std::flush;
                flush_waited = true;
            }
        }
        else {
            if (tty) std::cout << "\r";
            std::cout << "Waiting for output data to finish writing... (" <<
                creativity.storage().first->backend().pending() << " states pending)";
            if (tty) { std::cout << "   " << std::flush; need_endl = true; }
            else std::cout << std::endl;
        }
    }
    if (args.quiet and flush_waited) std::cout << "done" << std::endl;
    else if (need_endl) std::cout << std::endl;

    std::cout << "Compressing and saving final results file..." << std::flush;
    // Destroy the creativity object, which should cascade through to the storage object, causing it
    // to save its contents.
    cr_ptr.reset();
    std::cout << " done.\n\nSimulation saved to " << args_out << ".\n\n";
}
