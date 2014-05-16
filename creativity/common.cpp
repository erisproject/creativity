#include "creativity/common.hpp"
#include "creativity/Book.hpp"

namespace creativity {

eris::SharedMember<eris::Good::Continuous> MONEY;

std::vector<eris::SharedMember<Book>> NEW_BOOKS, AGE_ONE_BOOKS;

void NEW_BOOKS_Cleaner::intraFinish() {
    AGE_ONE_BOOKS.clear();
    NEW_BOOKS.swap(AGE_ONE_BOOKS);
}

}
