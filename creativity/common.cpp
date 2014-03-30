#include "creativity/common.hpp"
#include "creativity/BookMarket.hpp"

namespace creativity {

eris::SharedMember<eris::Good::Continuous> MONEY;

std::vector<eris::SharedMember<BookMarket>> NEW_BOOKS;

void NEW_BOOKS_Cleaner::intraFinish() { NEW_BOOKS.clear(); }

}
