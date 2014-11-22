#include "creativity/state/State.hpp"
#include "creativity/Reader.hpp"
#include "creativity/Book.hpp"

using namespace eris;

namespace creativity { namespace state {

State::State(const std::shared_ptr<Simulation> &sim) {
    auto lock = sim->runLock();

    t = sim->t();
    for (auto &r : sim->agents<Reader>()) {
        readers.emplace(
                std::piecewise_construct,
                std::tuple<eris_id_t>{r},
                std::tuple<Reader&>{r});
        if (boundary <= 0 or std::isnan(boundary)) boundary = r->wrapUpperBound()[0];
        if (dimensions == 0) dimensions = r->position().dimensions;
    }
    for (auto b : sim->goods<Book>()) {
        books.emplace(b->id(), b);
        if (boundary <= 0 or std::isnan(boundary)) boundary = b->wrapUpperBound()[0];
        if (dimensions == 0) dimensions = b->position().dimensions;
    }
}

}}
