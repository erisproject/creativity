#include "creativity/GUIShim.hpp"

namespace creativity {

GUIShim::GUIShim(Eris<Simulation> eris,
        std::function<void(GUI::Parameters)> &setup,
        std::function<void(unsigned int count)> &resume,
        std::function<void()> &stop) {
}

void GUIShim::start(int argc, char *argv[]) {
    gui_.start(argc, argv);
    gui_.join();
}

}
