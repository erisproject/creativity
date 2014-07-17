#include "creativity/state/ReaderState.hpp"
#include "creativity/Reader.hpp"

namespace creativity { namespace state {

ReaderState::ReaderState(const Reader &r) :
    id{r.id()},
    t{r.simulation()->t()},
    position{r.position()},
    u{r.u()},
    uLifetime{r.uLifetime()},
    belief{
        .profit = r.profitBelief(),
        .profit_extrap = r.profitExtrapBelief(),
        .demand = r.demandBelief(),
        .quality = r.qualityBelief(),
        .profit_stream = r.profitStreamBeliefs()
    }

{
    library.reserve(r.library().size());
    for (auto &bq : r.library()) {
        library.emplace(bq.first->id(), bq.second);
    }

    newBooks.reserve(r.newBooks().size());
    for (auto &b : r.newBooks()) {
        newBooks.emplace(b->id());
    }

    wrote.reserve(r.wrote().size());
    for (auto &b : r.wrote()) {
        wrote.emplace_back(b->id());
    }
}

}}
