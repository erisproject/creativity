#include "creativity/state/ReaderState.hpp"
#include "creativity/Reader.hpp"

namespace creativity { namespace state {

using namespace eris;

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
    cost_piracy{r.cost_piracy},
    income{r.income},
    profit{r.profitBelief()},
    profit_extrap{r.profitExtrapBelief()},
    demand{r.demandBelief()},
    quality{r.qualityBelief()},
    profit_stream{r.profitStreamBeliefs()}
{
    library.reserve(r.library().size());
    for (const auto &bq : r.library())
        library.emplace(bq.first->id(), bq.second);

    updateLibraryCounts(r.simulation()->t());

    friends.reserve(r.friends().size());
    for (const auto &m : r.friends()) friends.insert(m);

    new_books.reserve(r.newBooks().size());
    for (const auto &nb : r.newBooks()) new_books.insert(nb.first);

    for (const auto &b : r.wrote()) wrote.emplace_hint(wrote.end(), b->id());
}

void ReaderState::updateLibraryCounts(eris_time_t t) {
    library_purchased = 0;
    library_purchased_new = 0;
    library_pirated = 0;
    library_pirated_new = 0;
    for (auto &l : library) {
        if    (l.second.purchased()) { library_purchased++; if (l.second.acquired == t) library_purchased_new++; }
        else if (l.second.pirated()) { library_pirated++;   if (l.second.acquired == t) library_pirated_new++; }
    }
}

}}
