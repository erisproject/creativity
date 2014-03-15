#include "creativity/GUIShim.hpp"
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
    GUIShim gs(sim, setup, run, stop);

    try {
        gs.start(argc, argv);
    }
    catch (Glib::Error &e) {
        std::cerr << "Unable to start gui: " << e.what() << "\n";
        throw;
    }

}
