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
    for (const auto &bq : r.library())
        library.emplace(bq.first->id(), bq.second);

    copyIDs(r.friends(), friends);
    copyIDs(r.libraryPurchased(), library_purchased);
    copyIDs(r.libraryPirated(), library_pirated);
    copyIDs(r.newBooks(), new_books);
    copyIDs(r.newPurchased(), new_purchased);
    copyIDs(r.newPirated(), new_pirated);

    for (const auto &b : r.wrote()) wrote.emplace_hint(wrote.end(), b->id());
}

}}
