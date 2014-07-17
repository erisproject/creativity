#include "creativity/state/BookState.hpp"
#include "creativity/Book.hpp"
#include "creativity/Reader.hpp"

namespace creativity { namespace state {

BookState::BookState(const Book &b) : BookState(b, b.copies()) {}

BookState::BookState(const Book &b, unsigned long num_copies) :
    id{b.id()},
    t{b.simulation()->t()},
    author{b.author()->id()},
    position{b.position()},
    quality{b.quality()},
    market{b.hasMarket()},
    price{b.price()},
    revenue{b.currRevenue()},
    revenueLifetime{b.lifeRevenue()},
    sales{b.currSales()},
    salesLifetime{b.lifeSales()},
    copies{num_copies},
    age{b.age()},
    lifetime{b.marketPeriods()}
{}

}}
