#include "creativity/state/ReaderState.hpp"
#include "creativity/Reader.hpp"

namespace creativity { namespace state {

ReaderState::ReaderState(const unsigned int dimensions) :
    position{eris::Position::zero(dimensions)}
{}

ReaderState::ReaderState(const Reader &r) :
    id{r.id()},
    position{r.position()},
    u{r.u()},
    u_lifetime{r.uLifetime()},
    cost_fixed{r.cost_fixed},
    cost_unit{r.cost_unit},
    income{r.income},
    profit{r.profitBelief()},
    profit_extrap{r.profitExtrapBelief()},
    demand{r.demandBelief()},
    quality{r.qualityBelief()},
    profit_stream{r.profitStreamBeliefs()}
{
    library.reserve(r.library().size());
    for (auto &bq : r.library()) {
        library.emplace(bq.first->id(), bq.second);
    }

    new_books.reserve(r.newBooks().size());
    for (auto &b : r.newBooks()) {
        new_books.emplace(b->id());
    }

    wrote.reserve(r.wrote().size());
    for (auto &b : r.wrote()) {
        wrote.emplace_back(b->id());
    }
}

}}
