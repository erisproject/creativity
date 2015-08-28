#include "creativity/state/ReaderState.hpp"
#include <eris/Simulation.hpp>
#include <utility>
#include "creativity/Book.hpp"
#include "creativity/BookCopy.hpp"
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
    creation_shape{r.creation_shape},
    creation_scale{r.creation_scale},
    profit{r.profitBelief()},
    demand{r.demandBelief()},
    quality{r.qualityBelief()},
    profit_stream{r.profitStreamBeliefs()}
{
    if (r.profitExtrapBeliefDiffers()) profit_extrap = r.profitExtrapBelief();

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

const belief::Profit& ReaderState::profitExtrap() const {
    if (profit_extrap.n() > 0) return profit_extrap;
    return profit;
}

}}
