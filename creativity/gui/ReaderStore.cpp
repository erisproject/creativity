#include "creativity/gui/ReaderStore.hpp"
#include "creativity/gui/GUI.hpp"

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
    if (column == columns.id.index()) {
        Glib::Value<eris_id_t> v;
        v.init(v.value_type());
        v.set(r.id);
        value.init(v.gobj());
    }
    else if (column == columns.pos_x.index() or column == columns.pos_y.index() or column == columns.u.index() or column == columns.u_lifetime.index()
            or column == columns.cost_fixed.index() or column == columns.cost_unit.index() or column == columns.income.index()) {
        Glib::Value<double> v;
        v.init(v.value_type());
        v.set(  column == columns.pos_x.index() ? r.position[0] :
                column == columns.pos_y.index() ? r.position[1] :
                column == columns.u.index() ? r.u :
                column == columns.u_lifetime.index() ? r.u_lifetime :
                column == columns.cost_fixed.index() ? r.cost_fixed :
                column == columns.cost_unit.index() ? r.cost_unit :
                r.income
             );
        value.init(v.gobj());
    }
    else if (column == columns.pos_str.index()) {
        Glib::Value<std::string> v;
        v.init(v.value_type());
        v.set(GUI::pos_to_string(r.position));
        value.init(v.gobj());
    }
    else if (column == columns.books_owned.index() or column == columns.books_new.index() or column == columns.books_written.index()
            or column == columns.last_book_age.index() or column == columns.num_friends.index()
            or column == columns.books_pirated.index() or column == columns.books_purchased.index()
            or column == columns.books_new_pirated.index() or column == columns.books_new_purchased.index()) {
        Glib::Value<size_t> v;
        v.init(v.value_type());
        v.set(  column == columns.books_owned.index() ? r.library.size() :
                column == columns.books_purchased.index() ? r.library_purchased.size() :
                column == columns.books_pirated.index() ? r.library_pirated.size() :
                column == columns.books_new.index() ? r.new_books.size() :
                column == columns.books_new_purchased.index() ? r.new_purchased.size() :
                column == columns.books_new_pirated.index() ? r.new_pirated.size() :
                column == columns.books_written.index() ? r.wrote.size() :
                column == columns.num_friends.index() ? r.friends.size() :
                r.wrote.empty() ? state_->t : state_->books.at(*r.wrote.crbegin()).age
             );
        value.init(v.gobj());
    }
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
LESS_GREATER_A(books_purchased, library_purchased.size())
LESS_GREATER_A(books_pirated, library_pirated.size())
LESS_GREATER_A(books_new, new_books.size())
LESS_GREATER_A(books_new_purchased, new_purchased.size())
LESS_GREATER_A(books_new_pirated, new_pirated.size())
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
    return (a.wrote.empty() ? std::numeric_limits<unsigned long>::max() : state_->books.at(*a.wrote.crbegin()).age)
         < (b.wrote.empty() ? std::numeric_limits<unsigned long>::max() : state_->books.at(*b.wrote.crbegin()).age);
}
bool ReaderStore::greater_last_book_age(const ReaderState &a, const ReaderState &b) const {
    return (a.wrote.empty() ? std::numeric_limits<unsigned long>::max() : state_->books.at(*a.wrote.crbegin()).age)
         > (b.wrote.empty() ? std::numeric_limits<unsigned long>::max() : state_->books.at(*b.wrote.crbegin()).age);
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
