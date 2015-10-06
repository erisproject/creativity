#include "creativity/state/BookState.hpp"
#include "creativity/Book.hpp"
#include "creativity/Reader.hpp"
#include <eris/Position.hpp>
#include <eris/SharedMember.hpp>

namespace creativity { namespace state {

BookState::BookState(const unsigned int dimensions) :
    position{eris::Position::zero(dimensions)}
{}

BookState::BookState(const eris::SharedMember<Book> &b) : BookState((Book&) b) {}

BookState::BookState(const Book &b) :
    id{b.id()},
    author{b.author()->id()},
    position{b.position()},
    quality{b.qualityMean()},
    market_private{b.hasPrivateMarket()},
    price{b.price()},
    revenue{b.currRevenue()},
    revenue_lifetime{b.lifeRevenue()},
    prize{b.currPrize()},
    prize_lifetime{b.lifePrize()},
    sales{b.currSales()},
    sales_lifetime_private{b.lifeSalesPrivate()},
    sales_lifetime_public{b.lifeSalesPublic()},
    pirated{b.currPirated()},
    pirated_lifetime{b.lifePirated()},
    created{b.created()},
    lifetime_private{b.privateMarketPeriods()}
{}

}}
