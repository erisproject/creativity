#include "creativity/GUIShim.hpp"
#include "creativity/Book.hpp"
#include "creativity/Reader.hpp"

namespace creativity {

GUIShim::GUIShim(std::shared_ptr<Simulation> sim,
        std::function<void(GUI::Parameters)> &setup,
        std::function<void(unsigned int count)> &resume,
        std::function<void()> &stop)
    : sim_(sim)
{}

GUIShim::~GUIShim() {
    gui_.quit();
}

void GUIShim::start(int argc, char *argv[]) {
    gui_.start(argc, argv);
}

void GUIShim::sync() {
    gui_.clearPoints();
    gui_.clearCircles();
    for (auto &r : sim_->agents<Reader>()) {
        gui_.addPoint(r->position()[0], r->position()[1], GUI::PointType::X);
    }
    for (auto &b : sim_->goods<Book>()) {
        gui_.addPoint(b->position()[0], b->position()[1], GUI::PointType::SQUARE);
        gui_.addCircle(b->position()[0], b->position()[1], b->quality(), GUI::CircleType::A);
    }
    gui_.redraw();
}

}
