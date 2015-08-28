#include "creativity/state/State.hpp"
#include "creativity/Reader.hpp"
#include "creativity/Book.hpp"
#include <cmath>
#include <eris/Position.hpp>
#include <tuple>
#include <utility>

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
        if (dimensions == 0) dimensions = r->position().dimensions;
        if (boundary <= 0 or std::isnan(boundary)) boundary = r->wrapUpperBound()[0];
    }
    for (auto b : sim->goods<Book>()) {
        books.emplace(b->id(), b);
        if (dimensions == 0) dimensions = b->position().dimensions;
        if (boundary <= 0 or std::isnan(boundary)) boundary = b->wrapUpperBound()[0];
    }
}

}}
