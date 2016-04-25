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

int main(int argc, char *argv[]) {
    std::cerr << std::setprecision(16);
    std::cout << std::setprecision(16);

    auto cr_ptr = std::unique_ptr<Creativity>(new Creativity);
    auto &creativity = *cr_ptr;

    cmdargs::CLI args(creativity.set());
    args.parse(argc, argv);

    std::string args_out = args.output, tmp_out = args_out;
    if (args.xz) args_out = std::regex_replace(args_out, std::regex("(?:\\.[xX][zZ])?$"), ".xz");
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

        // We always write to an intermediate file, then move it into place at the end; by
        // default, that file is in the same directory as the output file, but the user might
        // want to put it somewhere else (e.g. if local storage is much faster than remote
        // storage, but the final result should be copied to the remote storage.)
        const auto &tmpdir = args.tmpdir;
        if (not tmpdir.empty()) {
            // First check and make sure that the parent of the requested output file exists
            if (not fs::is_directory(tmpdir))
                throw std::runtime_error("Error: --tmpdir `" + tmpdir + "' does not exist or is not a directory");

            dirname = tmpdir;
        }

        if (args.memory) {
            creativity.write<FileStorage>();
        }
        else {
            std::string tmpfile;
            int tries = 0;
            do {
                if (tries++ > 10) { std::cerr << "Error: unable to generate suitable tmpfile name (tried 10 times)\n"; exit(2); }
                tmpfile = dirname + "/" + random_filename(basename);
            } while (fs::exists(tmpfile));

            creativity.fileWrite(tmpfile);
            tmp_out = tmpfile;

            std::cout << "Writing simulation results to temp file: ";
            std::cout << tmp_out << "\n";
        }
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
                (creativity.publicSharing() ? u8"Public ✔; " : u8"Public ✘; ") <<
                "R=" << sim->countAgents<Reader>() << "; B=" << sim->countGoods<Book>() << "; Bnew=" <<
                std::setw(max_bnew_digits) << bnew << "] " << speed.str();
            if (tty) std::cout << std::flush;
            else std::cout << std::endl;
        }
    }

    if (not args.quiet) {
        if (tty) std::cout << std::endl;
        bool need_endl = false;
        while (not creativity.storage().first->flush_for(500)) {
            if (tty) std::cout << "\r";
            std::cout << "Waiting for output data to finish writing... (" <<
                creativity.storage().first->backend().pending() << " states pending)";
            if (tty) { std::cout << "   " << std::flush; need_endl = true; }
            else std::cout << std::endl;
        }
        if (need_endl) std::cout << std::endl;
    }
    else {
        creativity.storage().first->flush();
    }

    // Check again, in case something else created it while we were simulating:
    if (not overwrite) {
        if (fs::exists(args_out)) {
            std::cerr << "\nError: `" << args_out << "' already exists; specify a different file or add `--overwrite' option to overwrite\n";
            exit(1);
        }
    }

    std::string to_delete;
    if (args.xz) {
        std::cout << "Compressing..." << std::flush;
        try {
            dynamic_cast<FileStorage&>(creativity.storage().first->backend()).copyToXZ(args_out);
        }
        catch (const std::runtime_error &e) {
            std::cerr << "\nWriting to compressed file failed: " << e.what() << "\n\n";
            exit(1);
        }
        std::cout << " done." << std::endl;
        if (not args.memory) to_delete = tmp_out;
    }
    else if (args.memory) {
        std::cout << "Writing creativity data to " << args_out << "..." << std::flush;
        try {
            dynamic_cast<FileStorage&>(creativity.storage().first->backend()).copyTo(args_out);
        }
        catch (const std::runtime_error &e) {
            std::cerr << "\nWriting failed: " << e.what() << "\n\n";
            exit(1);
        }
        std::cout << " done." << std::endl;
    }
    else {
        // Destroy the Creativity object (so that the output file is closed)
        cr_ptr.reset();
        int error = std::rename(tmp_out.c_str(), args_out.c_str());
        if (error) { // Rename failed (perhaps across filesystems) so try a copy-and-delete

            std::cout << "Copying tmpfile to " << args_out << "..." << std::flush;
            {
                std::ifstream src; src.exceptions(src.failbit);
                std::ofstream dst; dst.exceptions(dst.failbit);
                try { src.open(tmp_out, std::ios::binary); }
                catch (std::ios_base::failure&) {
                    std::cerr << "\nUnable to read " << tmp_out << ": " << std::strerror(errno) << "\n";
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
            std::cout << " done.\nRemoving tmpfile..." << std::flush;

            // Copy succeeded, so delete the tmpfile
            to_delete = tmp_out;
        }
        tmp_out = args_out;
    }

    if (not to_delete.empty()) {
        int error = std::remove(to_delete.c_str());
        if (error) {
            // If we can't remove it, print a warning (but don't die, because we also copied it to the right place)
            std::cerr << "\nWarning: removing tmpfile `" << to_delete << "' failed: " << std::strerror(errno) << "\n";
        }
    }

    std::cout << "Simulation complete.  Results saved to " << args_out << "\n\n";
}
