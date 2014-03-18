#include "creativity/GUIShim.hpp"
#include "creativity/Reader.hpp"
#include "creativity/Book.hpp"
#include <eris/Eris.hpp>
#include <eris/Simulation.hpp>
#include <functional>
#include <iostream>

using namespace creativity;
using namespace std::placeholders;

class FakeHandlers {
    public:
        void setup(GUI::Parameters p) {
            std::cout << "FAKE SETUP HANDLER!\n";
        }
        void stop() {
            std::cout << "FAKE STOP HANDLER!\n";
        }
        void run(unsigned int rounds) {
            std::cout << "FAKE RUN HANDLER(" << rounds << ")!\n";
        }
};

int main(int argc, char *argv[]) {

    Eris<Simulation> sim;

    FakeHandlers f;

    std::function<void(GUI::Parameters)> setup = std::bind(&FakeHandlers::setup, f, _1);
    std::function<void()> stop = std::bind(&FakeHandlers::stop, f);
    std::function<void(unsigned int)> run = std::bind(&FakeHandlers::run, f, _1);
    GUIShim shim(sim, setup, run, stop);

    try {
        shim.start(argc, argv);
    }
    catch (Glib::Error &e) {
        std::cerr << "Unable to start gui: " << e.what() << "\n";
        throw;
    }

    for (long i = 0; i < 100000000000; i++) {
        if (i == 10000000000) {
            std::cout << "Adding agent at (-3,2), (-1,2), ..., (5,2)\n";
            sim->create<Reader>(Position{-3,2});
            sim->create<Reader>(Position{-1,2});
            sim->create<Reader>(Position{1,2});
            sim->create<Reader>(Position{3,2});
            sim->create<Reader>(Position{5,2});
            shim.sync();
        }
        if (i == 20000000000) {
            std::cout << "Adding book at (1,1,q=4) and (-3,0,q=2.5)\n";
            sim->create<Book>(Position{1,1},4);
            sim->create<Book>(Position{-3,0}, 2.5);
            shim.sync();
        }
    }
}
