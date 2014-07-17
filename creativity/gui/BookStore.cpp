#include "creativity/gui/BookStore.hpp"
#include "creativity/gui/GUI.hpp"
#include <gtkmm/treepath.h>
#include <algorithm>

using namespace eris;
using namespace std::placeholders;
using namespace creativity::state;

namespace creativity { namespace gui {
    
BookStore::BookStore(const State &state, eris_id_t author)
    : Glib::ObjectBase(typeid(BookStore)), MemberStore(state),
    author_{author}
{
    if (author) {
        auto &wrote = state_.readers.at(author).wrote;
        members_.reserve(wrote.size());
        for (auto &bid : wrote) {
            members_.emplace_back(state_.books.at(bid));
        }
    }
    else {
        members_.reserve(state_.books.size());
        for (auto &m : state_.books) {
            members_.emplace_back(m.second);
        }
    }
}

Glib::RefPtr<BookStore> BookStore::create(const State &state, eris_id_t author) {
    return Glib::RefPtr<BookStore>(new BookStore(state, author));
}

int BookStore::get_n_columns_vfunc() const {
    return columns.size();
}

GType BookStore::get_column_type_vfunc(int index) const {
    if ((unsigned int) index >= columns.size()) throw std::out_of_range("Invalid column index accessed");
    return columns.types()[index];
}

void BookStore::get_value_vfunc(const iterator &iter, int column, Glib::ValueBase &value) const {
    if (iter.get_stamp() != stamp_ or column > get_n_columns_vfunc()) return;

    auto &b = members_.at((size_t) iter.gobj()->user_data).get();
    if (column == columns.id.index() || column == columns.author.index()) {
        Glib::Value<eris_id_t> v;
        v.init(v.value_type());
        v.set(column == columns.id.index() ? b.id : b.author);
        value.init(v.gobj());
    }
    else if (column == columns.posX.index() or column == columns.posY.index() or column == columns.quality.index() or column == columns.price.index()
            or column == columns.revenue.index() or column == columns.revenueLifetime.index()) {
        Glib::Value<double> v;
        v.init(v.value_type());
        v.set(  column == columns.posX.index() ? b.position[0] :
                column == columns.posY.index() ? b.position[1] :
                column == columns.quality.index() ? b.quality :
                column == columns.price.index() ? b.price :
                column == columns.revenue.index() ? b.revenue :
                b.revenueLifetime
             );
        value.init(v.gobj());
    }
    else if (column == columns.market.index()) {
        Glib::Value<bool> v;
        v.init(v.value_type());
        v.set(b.market);
        value.init(v.gobj());
    }
    else if (column == columns.posstr.index()) {
        Glib::Value<std::string> v;
        v.init(v.value_type());
        v.set(GUI::pos_to_string(b.position));
        value.init(v.gobj());
    }
    else if (column == columns.age.index() or column == columns.sales.index() or column == columns.salesLifetime.index()
            or column == columns.lifetime.index() or column == columns.copies.index()) {
        Glib::Value<size_t> v;
        v.init(v.value_type());
        v.set(  column == columns.age.index() ? b.age :
                column == columns.sales.index() ? b.sales :
                column == columns.salesLifetime.index() ? b.salesLifetime :
                column == columns.lifetime.index() ? b.lifetime :
                b.copies
             );
        value.init(v.gobj());
    }
    else {
        throw std::out_of_range("Invalid column index accessed");
    }
}

void BookStore::set_sort_column_id_vfunc(int sort_column_id, Gtk::SortType order) {
    bool ascending = (order == Gtk::SORT_ASCENDING);
    std::function<bool(const BookState &a, const BookState &b)> compare;
#define ELSE_IF_COL(COL) \
    else if (sort_column_id == columns.COL.index()) compare = ascending ? less_##COL : greater_##COL
    if (sort_column_id == columns.id.index() || sort_column_id == DEFAULT_SORT_COLUMN_ID)
        compare = ascending ? less_id : greater_id;
    ELSE_IF_COL(author);
    ELSE_IF_COL(market);
    ELSE_IF_COL(posX);
    ELSE_IF_COL(posY);
    ELSE_IF_COL(posstr);
    ELSE_IF_COL(quality);
    ELSE_IF_COL(price);
    ELSE_IF_COL(revenue);
    ELSE_IF_COL(revenueLifetime);
    ELSE_IF_COL(age);
    ELSE_IF_COL(sales);
    ELSE_IF_COL(salesLifetime);
    ELSE_IF_COL(lifetime);
    ELSE_IF_COL(copies);
#undef ELSE_IF_COL

    sort_members(compare, sort_column_id, order);
}

#define LESS_GREATER(COL, ACCESS) \
bool BookStore::less_##COL(const BookState &a, const BookState &b) { \
    return a.ACCESS < b.ACCESS; \
} \
bool BookStore::greater_##COL(const BookState &a, const BookState &b) { \
    return a.ACCESS > b.ACCESS; \
}
LESS_GREATER(id, id)
LESS_GREATER(posX, position[0])
LESS_GREATER(posY, position[1])
LESS_GREATER(author, author)
LESS_GREATER(quality, quality)
LESS_GREATER(revenue, revenue)
LESS_GREATER(revenueLifetime, revenueLifetime)
LESS_GREATER(age, age)
LESS_GREATER(sales, sales)
LESS_GREATER(salesLifetime, salesLifetime)
LESS_GREATER(lifetime, lifetime)
LESS_GREATER(market, market)
LESS_GREATER(copies, copies)
#undef LESS_GREATER
// First x, then y for ties
bool BookStore::less_posstr(const BookState &a, const BookState &b) {
    auto ax = a.position[0], bx = b.position[0];
    return ax == bx ? a.position[1] < b.position[1] : ax < bx;
}
bool BookStore::greater_posstr(const BookState &a, const BookState &b) {
    auto ax = a.position[0], bx = b.position[0];
    return ax == bx ? a.position[1] > b.position[1] : ax > bx;
}
// For sorting purposes, marketless books are considered to have a price of negative infinity
// (rather than the model value of quiet_NaN) to sort them at one end.
bool BookStore::less_price(const BookState &a, const BookState &b) {
    return (a.market ? a.price : -std::numeric_limits<double>::infinity())
         < (b.market ? b.price : -std::numeric_limits<double>::infinity());
}
bool BookStore::greater_price(const BookState &a, const BookState &b) {
    return (a.market ? a.price : -std::numeric_limits<double>::infinity())
         > (b.market ? b.price : -std::numeric_limits<double>::infinity());
}

void BookStore::appendColumnsTo(Gtk::TreeView &v) const {
    appendCol(v, "ID", columns.id, 50);
    appendCol(v, "Position", columns.posstr, 100);
    if (not author_) appendCol(v, "Author", columns.author, 80);
    appendCol(v, "Age", columns.age, 65);
    appendCol(v, "Quality", columns.quality, 85);
    appendCol(v, "Mkt?", columns.market, 65);
    appendCol(v, "Price", columns.price, 75);
    appendCol(v, "Rev.", columns.revenueLifetime, 75);
    appendCol(v, "Sales", columns.salesLifetime, 75);
    appendCol(v, "Copies", columns.copies, 80);
    appendCol(v, "Life", columns.lifetime, 65);
}

}}
