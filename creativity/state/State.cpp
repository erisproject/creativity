#include "creativity/state/State.hpp"
#include "creativity/Reader.hpp"
#include "creativity/Book.hpp"

using namespace eris;

namespace creativity { namespace state {

State::State(const std::shared_ptr<Simulation> &sim) {
    auto lock = sim->runLock();

    t = sim->t();
    // We can get the number of copies of books a bit more efficiently than calling book->copies()
    // for every book (the latter requires a loop through all readers for all books, thus requiring
    // R*B iterations (where R=# readers, B=# books); this way requires R*b iterations, where b is
    // the average number of books owned per reader, which will typically be much smaller than B.
    std::unordered_map<eris_id_t, unsigned long> copies;
    for (auto r : sim->agents<Reader>()) {
        readers.emplace(
                std::piecewise_construct,
                std::tuple<eris_id_t>{r},
                std::tuple<Reader&>{r});
        for (auto b : r->library()) {
            copies[b.first]++;
        }
        if (boundary == 0) boundary = r->wrapUpperBound()[0];
        if (dimensions == 0) dimensions = r->position().dimensions;
    }
    for (auto b : sim->goods<Book>()) {
        books.emplace(
                std::piecewise_construct,
                std::tuple<eris_id_t>{b},
                std::tuple<Book&, unsigned long>{b, copies[b]});
        if (boundary == 0) boundary = b->wrapUpperBound()[0];
        if (dimensions == 0) dimensions = b->position().dimensions;
    }
}

}}
