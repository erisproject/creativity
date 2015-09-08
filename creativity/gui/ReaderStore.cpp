#include "creativity/gui/ReaderStore.hpp"
#include "creativity/gui/GUI.hpp"
#include <eris/Position.hpp>
#include <functional>
#include <iterator>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <glibmm/objectbase.h>
#include <glibmm/value.h> // IWYU pragma: keep
#include <glibmm/refptr.h>
#include <gtkmm/enums.h>
#include <gtkmm/treesortable.h>
#include <gtkmm/treemodel.h>
#include <gtkmm/treemodelcolumn.h>
#include <gtkmm/treeview.h> // IWYU pragma: keep
#include <gtkmm/treepath.h>


using namespace eris;
using namespace creativity::state;
using namespace std::placeholders;

namespace creativity { namespace gui {

    
ReaderStore::ReaderStore(std::shared_ptr<const State> &&state) : Glib::ObjectBase(typeid(ReaderStore)), MemberStore(std::move(state)) {
    members_.reserve(state_->readers.size());
    for (auto &m : state_->readers) {
        members_.emplace_back(m.second);
    }
}

Glib::RefPtr<ReaderStore> ReaderStore::create(std::shared_ptr<const State> s) {
    return Glib::RefPtr<ReaderStore>(new ReaderStore(std::move(s)));
}

int ReaderStore::get_n_columns_vfunc() const {
    return columns.size();
}

GType ReaderStore::get_column_type_vfunc(int index) const {
    if ((unsigned int) index >= columns.size()) throw std::out_of_range("Invalid column index accessed");
    return columns.types()[index];
}

void ReaderStore::get_value_vfunc(const iterator &iter, int column, Glib::ValueBase &value) const {
    if (iter.get_stamp() != stamp_ or column > get_n_columns_vfunc()) return;

    auto &r = members_.at((size_t) iter.gobj()->user_data).get();

// If statement for field with arbitrary value:
#define IFCOL_V(F, V) if (column == columns.F.index()) copy_value_(value, V)
// If statement for field with member of same name as the field:
#define IFCOL(F) IFCOL_V(F, r.F)
    IFCOL(id);
    else IFCOL_V(pos_x, r.position[0]);
    else IFCOL_V(pos_y, r.position[1]);
    else IFCOL(u);
    else IFCOL(u_lifetime);
    else IFCOL_V(pos_str, GUI::pos_to_string(r.position));
    else IFCOL_V(books_owned, r.library.size());
    else IFCOL_V(books_purchased, r.library_purchased);
    else IFCOL_V(books_pirated, r.library_pirated);
    else IFCOL_V(books_new, r.new_books.size());
    else IFCOL_V(books_new_purchased, r.library_purchased_new);
    else IFCOL_V(books_new_pirated, r.library_pirated_new);
    else IFCOL_V(books_written, r.wrote.size());
    else IFCOL_V(num_friends, r.friends.size());
    else IFCOL_V(last_book_age, r.wrote.empty() ? state_->t : state_->t - state_->books.at(*r.wrote.crbegin()).created);
#undef IFCOL
#undef IFCOL_V
    else {
        throw std::out_of_range("Invalid column index accessed");
    }
}

void ReaderStore::set_sort_column_id_vfunc(int sort_column_id, Gtk::SortType order) {
    bool ascending = (order == Gtk::SORT_ASCENDING);
    std::function<bool(const ReaderState &a, const ReaderState &b)> compare;
#define ELSE_IF_COL(COL) \
    else if (sort_column_id == columns.COL.index()) compare = ascending ? less_##COL : greater_##COL
    if (sort_column_id == columns.id.index() || sort_column_id == DEFAULT_SORT_COLUMN_ID)
        compare = ascending ? less_id : greater_id;
    ELSE_IF_COL(pos_x);
    ELSE_IF_COL(pos_y);
    ELSE_IF_COL(pos_str);
    ELSE_IF_COL(u);
    ELSE_IF_COL(u_lifetime);
    ELSE_IF_COL(num_friends);
    ELSE_IF_COL(books_owned);
    ELSE_IF_COL(books_purchased);
    ELSE_IF_COL(books_pirated);
    ELSE_IF_COL(books_new);
    ELSE_IF_COL(books_new_purchased);
    ELSE_IF_COL(books_new_pirated);
    ELSE_IF_COL(books_written);
#undef ELSE_IF_COL
    else if (sort_column_id == columns.last_book_age.index())
        // Can't handle like the above: less/greater_copies() is not static
        compare = std::bind(ascending ? &ReaderStore::less_last_book_age : &ReaderStore::greater_last_book_age, this, _1, _2);

    sort_members(compare, sort_column_id, order);
}


#define LESS_GREATER_A(COL, ACCESS) \
bool ReaderStore::less_##COL   (const ReaderState &a, const ReaderState &b) { return a.ACCESS < b.ACCESS; } \
bool ReaderStore::greater_##COL(const ReaderState &a, const ReaderState &b) { return a.ACCESS > b.ACCESS; }
#define LESS_GREATER(FIELD) LESS_GREATER_A(FIELD, FIELD)
LESS_GREATER(id)
LESS_GREATER_A(pos_x, position[0])
LESS_GREATER_A(pos_y, position[1])
LESS_GREATER(u)
LESS_GREATER(u_lifetime)
LESS_GREATER_A(num_friends, friends.size())
LESS_GREATER_A(books_owned, library.size())
LESS_GREATER_A(books_purchased, library_purchased)
LESS_GREATER_A(books_pirated, library_pirated)
LESS_GREATER_A(books_new, new_books.size())
LESS_GREATER_A(books_new_purchased, library_purchased_new);
LESS_GREATER_A(books_new_pirated, library_pirated_new);
LESS_GREATER_A(books_written, wrote.size())
#undef LESS_GREATER
#undef LESS_GREATER_A

// First x, then y for ties
bool ReaderStore::less_pos_str(const ReaderState &a, const ReaderState &b) {
    auto ax = a.position[0], bx = b.position[0];
    return ax == bx ? a.position[1] < b.position[1] : ax < bx;
}
bool ReaderStore::greater_pos_str(const ReaderState &a, const ReaderState &b) {
    auto ax = a.position[0], bx = b.position[0];
    return ax == bx ? a.position[1] > b.position[1] : ax > bx;
}
bool ReaderStore::less_last_book_age(const ReaderState &a, const ReaderState &b) const {
    return (a.wrote.empty() ? 0 : state_->books.at(*a.wrote.crbegin()).created)
         > (b.wrote.empty() ? 0 : state_->books.at(*b.wrote.crbegin()).created);
}
bool ReaderStore::greater_last_book_age(const ReaderState &a, const ReaderState &b) const {
    return (a.wrote.empty() ? 0 : state_->books.at(*a.wrote.crbegin()).created)
         < (b.wrote.empty() ? 0 : state_->books.at(*b.wrote.crbegin()).created);
}

void ReaderStore::appendColumnsTo(Gtk::TreeView &v) const {
    appendCol(v, "ID", columns.id, 50);
    appendCol(v, "Position", columns.pos_str, 110);
    appendCol(v, "Utility", columns.u, 90);
    appendCol(v, "Life Util.", columns.u_lifetime, 90);
    appendCol(v, "# Friends", columns.num_friends, 80);
    appendCol(v, "Books", columns.books_owned, 80);
    appendCol(v, "# New", columns.books_new, 80);
    appendCol(v, "# Bought", columns.books_purchased, 100);
    appendCol(v, "# Pirated", columns.books_pirated, 100);
    appendCol(v, "# Written", columns.books_written, 100);
    appendCol(v, "Last wrote", columns.last_book_age, 100);
}

}}
