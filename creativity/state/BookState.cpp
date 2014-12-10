#include "creativity/state/BookState.hpp"
#include "creativity/Book.hpp"
#include "creativity/Reader.hpp"

namespace creativity { namespace state {

BookState::BookState(const unsigned int dimensions) :
    position{eris::Position::zero(dimensions)}
{}

BookState::BookState(const eris::SharedMember<Book> &b) : BookState((Book&) b) {}

BookState::BookState(const Book &b) :
    id{b.id()},
    author{b.author()->id()},
    position{b.position()},
    quality{b.quality()},
    price{b.price()},
    revenue{b.currRevenue()},
    revenue_lifetime{b.lifeRevenue()},
    sales{b.currSales()},
    sales_lifetime{b.lifeSales()},
    pirated{b.currPirated()},
    pirated_lifetime{b.lifePirated()},
    created{b.created()},
    lifetime{b.marketPeriods()}
{}

}}
