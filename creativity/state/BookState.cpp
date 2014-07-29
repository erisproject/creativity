#include "creativity/state/BookState.hpp"
#include "creativity/Book.hpp"
#include "creativity/Reader.hpp"

namespace creativity { namespace state {

BookState::BookState(const unsigned int dimensions) :
    position{eris::Position::zero(dimensions)}
{}

BookState::BookState(const Book &b) : BookState(b, b.copies()) {}

BookState::BookState(const Book &b, unsigned long num_copies) :
    id{b.id()},
    author{b.author()->id()},
    position{b.position()},
    quality{b.quality()},
    market{b.hasMarket()},
    price{b.price()},
    revenue{b.currRevenue()},
    revenue_lifetime{b.lifeRevenue()},
    sales{b.currSales()},
    sales_lifetime{b.lifeSales()},
    copies{num_copies},
    age{b.age()},
    lifetime{b.marketPeriods()}
{}

}}
