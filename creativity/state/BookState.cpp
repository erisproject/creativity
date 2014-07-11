#include "creativity/state/BookState.hpp"
#include "creativity/Book.hpp"
#include "creativity/Reader.hpp"

namespace creativity { namespace state {

BookState::BookState(const Book &b) :
    id{b.id()},
    t{b.simulation()->t()},
    author{b.author()->id()},
    position{b.position()},
    quality{b.quality()},
    market{b.hasMarket()},
    price{b.hasMarket() ? b.market()->price() : std::numeric_limits<double>::quiet_NaN()},
    revenue{b.currRevenue()},
    revenueLifetime{b.lifeRevenue()},
    sales{b.currSales()},
    salesLifetime{b.lifeSales()},
    copies{b.copies()},
    age{b.age()},
    lifetime{b.marketPeriods()}
{}

}}
