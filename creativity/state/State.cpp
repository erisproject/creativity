#include "creativity/state/State.hpp"
#include "creativity/Reader.hpp"
#include "creativity/Book.hpp"

using namespace eris;

namespace creativity { namespace state {

State::State(const std::shared_ptr<Simulation> &sim) {
    auto lock = sim->runLock();

    // Bypass the const modifiers here in the constructor:
    const_cast<unsigned long&>(t) = sim->t();
    auto &rdrs = const_cast<std::unordered_map<eris_id_t, ReaderState>&>(readers);
    auto &bks = const_cast<std::unordered_map<eris_id_t, BookState>&>(books);

    // We can get the number of copies of books a bit more efficiently than calling book->copies()
    // for every book (the latter requires a loop through all readers for all books, thus requiring
    // R*B iterations (where R=# readers, B=# books); this way requires R*b iterations, where b is
    // the average number of books owned per reader, which will typically be much smaller than B.
    std::unordered_map<eris_id_t, unsigned long> copies;
    for (auto r : sim->agents<Reader>()) {
        rdrs.emplace(
                std::piecewise_construct,
                std::tuple<eris_id_t>{r},
                std::tuple<decltype(r)>{r});
        for (auto b : r->library()) {
            copies[b.first]++;
        }
    }
    for (auto b : sim->goods<Book>()) {
        bks.emplace(
                std::piecewise_construct,
                std::tuple<eris_id_t>{b},
                std::tuple<decltype(b), unsigned long>{b, copies[b]});
    }
}

}}
