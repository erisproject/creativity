#pragma once
#include <eris/SharedMember.hpp>
#include <eris/Good.hpp>
#include <eris/Optimize.hpp>

/// \file creativity/common.hpp Defines various global objects and variables, such as the money good.

namespace creativity {

/// The money good.  Set when created.
extern eris::SharedMember<eris::Good::Continuous> MONEY;

/** The reader boundary (from 0, in each dimension).  If readers attempt to move off the boundary,
 * they "wrap", coming in the opposite side(s).
 */
constexpr double BOUNDARY = 10;

class BookMarket; // Forward declaration
/** List of brand-new BookMarkets created this period.  Cleared at the end of each
 * intra-optimization stage, built up when books are created (in Book::interApply())
 */
extern std::vector<eris::SharedMember<BookMarket>> NEW_BOOKS;

/** Simple class that clears NEW_BOOKS at the end of every period.  There should be one and only one
 * instance of this class in the simulation.
 */
class NEW_BOOKS_Cleaner : public eris::Member, public virtual eris::intraopt::Finish {
    public:
        /// Clears NEW_BOOKS at the end of a period.
        void intraFinish() override;
};

}
