#include "creativity/gui/BookStore.hpp"
#include "creativity/gui/GUI.hpp"
#include <gtkmm/treepath.h>
#include <algorithm>

using namespace eris;
using namespace std::placeholders;
using namespace creativity::state;

namespace creativity { namespace gui {

BookStore::BookStore(std::shared_ptr<const State> &&state, eris_id_t author, std::unique_ptr<ColRec> &&cols)
    : Glib::ObjectBase(typeid(BookStore)), MemberStore(std::move(state)),
    columns(std::move(cols)), author_(std::move(author))
{
    initializeBooks();
}

BookStore::BookStore(std::shared_ptr<const State> &&state, eris_id_t author)
    : Glib::ObjectBase(typeid(BookStore)), MemberStore(std::move(state)), author_(std::move(author))
{
    initializeBooks();
}

BookStore::BookStore(std::shared_ptr<const State> &&state, std::unique_ptr<ColRec> &&cols)
    : Glib::ObjectBase(typeid(BookStore)), MemberStore(std::move(state)), columns(std::move(cols))
{
    // Don't initializeBooks()--the subclass is going to take care of it
}


Glib::RefPtr<BookStore> BookStore::create(std::shared_ptr<const State> state, eris_id_t author) {
    return Glib::RefPtr<BookStore>(new BookStore(std::move(state), author));
}

void BookStore::initializeBooks() {
    if (author_) {
        auto &wrote = state_->readers.at(author_).wrote;
        members_.reserve(wrote.size());
        for (auto &bid : wrote) {
            members_.emplace_back(state_->books.at(bid));
        }
    }
    else {
        members_.reserve(state_->books.size());
        for (auto &m : state_->books) {
            members_.emplace_back(m.second);
        }
    }
}

int BookStore::get_n_columns_vfunc() const {
    return columns->size();
}

GType BookStore::get_column_type_vfunc(int index) const {
    if ((unsigned int) index >= columns->size()) throw std::out_of_range("Invalid column index accessed");
    return columns->types()[index];
}

void BookStore::get_value_vfunc(const iterator &iter, int column, Glib::ValueBase &value) const {
    if (iter.get_stamp() != stamp_ or column > get_n_columns_vfunc()) return;

    auto &b = members_.at((size_t) iter.gobj()->user_data).get();
// If statement for field with arbitrary value:
#define IFCOL_V(F, V) if (column == columns->F.index()) copy_value_(value, V)
// If statement for field with member of same name as the field:
#define IFCOL(F) IFCOL_V(F, b.F)
// If statement for field with member method of same name as the field:
#define IFCOL_M(F) IFCOL_V(F, b.F())
    IFCOL(id);
    else IFCOL(author);
    else IFCOL_V(pos_x, b.position[0]);
    else IFCOL_V(pos_y, b.position[1]);
    else IFCOL(quality);
    else IFCOL(price);
    else IFCOL(revenue);
    else IFCOL(revenue_lifetime);
    else IFCOL_V(pos_str, GUI::pos_to_string(b.position));
    else IFCOL_V(market_str, b.market_private ? "Priv." : b.market_any() ? "Pub." : "No");
    else IFCOL_M(market_any);
    else IFCOL(market_private);
    else IFCOL_M(market_public);
    else IFCOL_V(age, state_->t - b.created);
    else IFCOL(created);
    else IFCOL(sales);
    else IFCOL(sales_lifetime_private);
    else IFCOL(sales_lifetime_public);
    else IFCOL_M(sales_lifetime);
    else IFCOL(pirated);
    else IFCOL(pirated_lifetime);
    else IFCOL_M(copies);
    else IFCOL(lifetime_private);
    else throw std::out_of_range("Invalid column index accessed");
#undef IFCOL_M
#undef IFCOL
#undef IFCOL_V
}

void BookStore::set_sort_column_id_vfunc(int sort_column_id, Gtk::SortType order) {
    bool ascending = (order == Gtk::SORT_ASCENDING);
    std::function<bool(const BookState &a, const BookState &b)> compare;
    if (sort_column_id == columns->id.index() || sort_column_id == DEFAULT_SORT_COLUMN_ID)
        compare = ascending ? less_id : greater_id;
#define ELSE_IF_COL(COL) \
    else if (sort_column_id == columns->COL.index()) compare = ascending ? less_##COL : greater_##COL
    ELSE_IF_COL(author);
    ELSE_IF_COL(pos_x);
    ELSE_IF_COL(pos_y);
    ELSE_IF_COL(quality);
    ELSE_IF_COL(price);
    ELSE_IF_COL(revenue);
    ELSE_IF_COL(revenue_lifetime);
    ELSE_IF_COL(pos_str);
    ELSE_IF_COL(market_private);
    ELSE_IF_COL(market_public);
    ELSE_IF_COL(market_any);
    ELSE_IF_COL(age);
    ELSE_IF_COL(created);
    ELSE_IF_COL(sales);
    ELSE_IF_COL(sales_lifetime_private);
    ELSE_IF_COL(sales_lifetime_public);
    ELSE_IF_COL(sales_lifetime);
    ELSE_IF_COL(pirated);
    ELSE_IF_COL(pirated_lifetime);
    ELSE_IF_COL(copies);
    ELSE_IF_COL(lifetime_private);
#undef ELSE_IF_COL

    sort_members(compare, sort_column_id, order);
}

#define LESS_GREATER_A(COL, ACCESS) \
bool BookStore::less_##COL   (const BookState &a, const BookState &b) { return a.ACCESS < b.ACCESS; } \
bool BookStore::greater_##COL(const BookState &a, const BookState &b) { return a.ACCESS > b.ACCESS; }
#define LESS_GREATER(FIELD) LESS_GREATER_A(FIELD, FIELD)
#define LESS_GREATER_M(FIELD) LESS_GREATER_A(FIELD, FIELD())
LESS_GREATER(id)
LESS_GREATER(author)
LESS_GREATER_A(pos_x, position[0])
LESS_GREATER_A(pos_y, position[1])
LESS_GREATER(quality)
// price handled below
LESS_GREATER(revenue)
LESS_GREATER(revenue_lifetime)
// pos_str handled below
// market_str handled below
LESS_GREATER(market_private)
LESS_GREATER_M(market_public)
LESS_GREATER_M(market_any)
// age handled below
LESS_GREATER(created)
LESS_GREATER(sales)
LESS_GREATER(sales_lifetime_private)
LESS_GREATER(sales_lifetime_public)
LESS_GREATER_M(sales_lifetime)
LESS_GREATER(pirated)
LESS_GREATER(pirated_lifetime)
LESS_GREATER_A(copies, copies_lifetime())
LESS_GREATER(lifetime_private)
#undef LESS_GREATER
#undef LESS_GREATER_A
// First x, then y for ties
bool BookStore::less_pos_str(const BookState &a, const BookState &b) {
    auto ax = a.position[0], bx = b.position[0];
    return ax == bx ? a.position[1] < b.position[1] : ax < bx;
}
bool BookStore::greater_pos_str(const BookState &a, const BookState &b) {
    auto ax = a.position[0], bx = b.position[0];
    return ax == bx ? a.position[1] > b.position[1] : ax > bx;
}
// Sorting by market: No < Pub < Priv
bool BookStore::less_market_str(const BookState &a, const BookState &b) {
    return
        (a.market_private ? 2 : a.market_any() ? 1 : 0)
        <
        (b.market_private ? 2 : b.market_any() ? 1 : 0);
}
bool BookStore::greater_market_str(const BookState &a, const BookState &b) {
    return
        (a.market_private ? 2 : a.market_any() ? 1 : 0)
        >
        (b.market_private ? 2 : b.market_any() ? 1 : 0);
}

// For sorting purposes, marketless books are considered to have a price of negative infinity
// (rather than the model value of quiet_NaN) to sort them at one end.
bool BookStore::less_price(const BookState &a, const BookState &b) {
    return (a.market_any() ? a.price : -std::numeric_limits<double>::infinity())
         < (b.market_any() ? b.price : -std::numeric_limits<double>::infinity());
}
bool BookStore::greater_price(const BookState &a, const BookState &b) {
    return (a.market_any() ? a.price : -std::numeric_limits<double>::infinity())
         > (b.market_any() ? b.price : -std::numeric_limits<double>::infinity());
}
// Age sorts in the opposite order from created
bool BookStore::less_age   (const BookState &a, const BookState &b) { return a.created > b.created; }
bool BookStore::greater_age(const BookState &a, const BookState &b) { return a.created < b.created; }

void BookStore::appendColumnsTo(Gtk::TreeView &v) const {
    appendCol(v, "ID", columns->id, 50);
    appendCol(v, "Position", columns->pos_str, 110);
    if (not author_) appendCol(v, "Author", columns->author, 80);
    appendCol(v, "Created", columns->created, 85);
    appendCol(v, "Life", columns->lifetime_private, 65);
    appendCol(v, "Mkt", columns->market_str, 65);
    appendCol(v, "Quality", columns->quality, 85);
    appendCol(v, "Price", columns->price, 75);
    appendCol(v, "Rev.", columns->revenue_lifetime, 75);
    appendCol(v, "Sales", columns->sales_lifetime, 75);
    appendCol(v, "Pirated", columns->pirated_lifetime, 85);
    appendCol(v, "Copies", columns->copies, 80);
}

}}
