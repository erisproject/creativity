#include "creativity/Creativity.hpp"
#include "creativity/Reader.hpp"
#include "creativity/Book.hpp"
#include "creativity/state/Storage.hpp"
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
#include <cstdio>
#include <boost/filesystem/path.hpp>
#include <sys/stat.h>

using namespace creativity;
using namespace creativity::state;
using namespace eris;


bool exists_dir(const std::string &p) {
    struct stat buffer;
    return (stat(p.c_str(), &buffer) == 0 and S_ISDIR(buffer.st_mode));
}

bool exists(const std::string &p) {
    struct stat buffer;
    return (stat(p.c_str(), &buffer) == 0);
}

std::string random_filename(const std::string &basename) {
    std::random_device rd;
    std::uniform_int_distribution<unsigned char> rand_index(0, 35);
    constexpr char map[37] = "0123456789abcdefghijklmnopqrstuvwxyz";

    std::ostringstream buf;
    buf << basename << ".partial.";
    for (int i = 0; i < 7; i++) {
        buf << map[rand_index(rd)];
    }
    return buf.str();
}

int main(int argc, char *argv[]) {
    std::cerr << std::setprecision(16);
    std::cout << std::setprecision(16);
    auto creativity = Creativity::create();

    cmdargs::CLI cmd(creativity->set());
    cmd.parse(argc, argv);

    std::string results_out = cmd.output;
    std::string args_out = results_out;
    bool overwrite = cmd.overwrite;
    bool need_copy = false; // Whether we need to copy from args_out to results_out (if using --tmpdir)

    // Filename output

    try {
        // If the user didn't specify --overwrite, make sure the file doesn't exist.  (It might
        // get created before we finish the simulation, but at least we can abort now in a
        // typical case).
        if (not overwrite) {
            if (exists(args_out)) {
                std::cerr << "Error: `" << args_out << "' already exists; specify a different file or add `--overwrite' option to overwrite\n";
                exit(1);
            }
        }


        boost::filesystem::path outpath(args_out);
        std::string basename = outpath.filename().string();
        std::string dirname = outpath.parent_path().string();
        if (dirname.empty()) dirname = ".";

        if (basename.empty()) throw std::runtime_error("Invalid output filename (" + args_out + ")");

        // Make sure the parent of the output file exists
        if (not exists_dir(dirname))
            throw std::runtime_error("Directory `" + dirname + "' does not exist or is not a directory");

        // We always write to an intermediate file, then move it into place at the end; by
        // default, that file is in the same directory as the output file, but the user might
        // want to put it somewhere else (e.g. if local storage is much faster than remote
        // storage, but the final result should be copied to the remote storage.)
        const auto &tmpdir = cmd.tmpdir;
        if (not tmpdir.empty()) {
            // First check and make sure that the parent of the requested output file exists
            if (not exists_dir(tmpdir))
                throw std::runtime_error("Error: --tmpdir `" + tmpdir + "' does not exist or is not a directory");

            dirname = tmpdir;
        }

        std::string tmpfile;
        int tries = 0;
        do {
            if (tries++ > 10) { std::cerr << "Error: unable to generate suitable tmpfile name (tried 10 times)\n"; exit(2); }
            tmpfile = dirname + "/" + random_filename(basename);
        } while (exists(tmpfile));

        creativity->fileWrite(tmpfile);
        results_out = tmpfile;
        need_copy = true;
    }
    catch (std::exception &e) {
        std::cerr << "Unable to write to file: " << e.what() << "\n";
        exit(1);
    }
    std::cout << "Writing simulation results to ";
    if (need_copy) std::cout << "temp file: ";
    std::cout << results_out << "\n";

    creativity->setup();
    auto sim = creativity->sim;

    // Copy the initial state into the storage object
    creativity->storage().first->emplace_back(sim);

    if (cmd.threads > 0) Eigen::initParallel();
    sim->maxThreads(cmd.threads);

    constexpr size_t avg_times = 5;
    std::chrono::time_point<std::chrono::high_resolution_clock> started = std::chrono::high_resolution_clock::now();
    std::queue<std::chrono::time_point<std::chrono::high_resolution_clock>> times;
    times.push(started);
    unsigned int max_bnew_digits = 0;
    auto count_bnew = [](const Book &b) -> bool { return b.age() == 0; };

    unsigned int periods = cmd.periods;

    while (sim->t() < periods) {
        creativity->run();

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


        std::cout << "\rRunning simulation [t=" << sim->t() << "; " <<
            (sim->t() >= creativity->parameters.piracy_begins ? u8"P✔" : u8"P✘") <<
            "; R=" << sim->countAgents<Reader>() << "; B=" << sim->countGoods<Book>() << "; Bnew=" <<
            std::setw(max_bnew_digits) << bnew << "] " << speed.str() << std::flush;
    }

    std::cout << "\nSimulation done.";
    unsigned int pending = creativity->storage().first->backend().pending();
    if (pending > 0) {
        std::string waiting = "Waiting for output data to finish writing... ";
        std::cout << "\n" << waiting;
        while (pending > 0) {
            std::cout << "\r" << waiting << "(" << pending << " states pending)  " << std::flush;
            std::this_thread::sleep_for(std::chrono::milliseconds(250));
            pending = creativity->storage().first->backend().pending();
        }
        std::cout << "\r" << waiting << "done.               " << std::flush;
    }

    std::cout << "\n";
    creativity->storage().first->flush();

    // Destroy the simulation (which also closes the storage)
    creativity.reset();

    if (need_copy) {
        std::cout << "Moving tmpfile to " << args_out << "..." << std::flush;
        int error = std::rename(results_out.c_str(), args_out.c_str());
        if (error) { // Rename failed (perhaps across filesystems) so do a copy-and-delete
            if (not overwrite) {
                if (exists(args_out)) {
                    std::cerr << "\nError: `" << args_out << "' already exists; specify a different file or add `--overwrite' option to overwrite\n";
                    exit(1);
                }
            }

            {
                std::ifstream src; src.exceptions(src.failbit);
                std::ofstream dst; dst.exceptions(dst.failbit);
                try { src.open(results_out, std::ios::binary); }
                catch (std::ios_base::failure&) {
                    std::cerr << "\nUnable to read " << results_out << ": " << std::strerror(errno) << "\n";
                    exit(1);
                }

                try { dst.open(args_out, std::ios::binary | std::ios::trunc); }
                catch (std::ios_base::failure&) {
                    std::cerr << "\nUnable to write to " << args_out << ": " << std::strerror(errno) << "\n";
                    exit(1);
                }

                try { dst << src.rdbuf(); }
                catch (std::ios_base::failure&) {
                    std::cerr << "\nUnable to copy file contents: " << std::strerror(errno) << "\n";
                    exit(1);
                }
            }

            std::cout << " done.\n" << std::flush;

            // Copy succeeded, so delete the tmpfile
            error = std::remove(results_out.c_str());
            if (error) {
                // If we can't remove it, print a warning (but don't die, because we also copied it to the right place)
                std::cerr << "\nWarning: removing tmpfile `" << results_out << "' failed: " << std::strerror(errno) << "\n";
            }
        }
        else {
            std::cout << " done.\n" << std::flush;
        }
        results_out = args_out;
    }

    std::cout << "\nSimulation complete.  Results saved to " << results_out << "\n\n";
}
