#include "creativity/gui/ReaderStore.hpp"
#include "creativity/gui/GUI.hpp"

using namespace eris;
using namespace creativity::state;
using namespace std::placeholders;

namespace creativity { namespace gui {

    
ReaderStore::ReaderStore(const State &state) : Glib::ObjectBase(typeid(ReaderStore)), MemberStore(state) {
    members_.reserve(state_.readers.size());
    for (auto &m : state_.readers) {
        members_.emplace_back(m.second);
    }
}

Glib::RefPtr<ReaderStore> ReaderStore::create(const State &s) {
    return Glib::RefPtr<ReaderStore>(new ReaderStore(s));
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
    else if (column == columns.posX.index() or column == columns.posY.index() or column == columns.u.index() or column == columns.uLifetime.index()) {
        Glib::Value<double> v;
        v.init(v.value_type());
        v.set(  column == columns.posX.index() ? r.position[0] :
                column == columns.posY.index() ? r.position[1] :
                column == columns.u.index() ? r.u :
                r.uLifetime
             );
        value.init(v.gobj());
    }
    else if (column == columns.posstr.index()) {
        Glib::Value<std::string> v;
        v.init(v.value_type());
        v.set(GUI::pos_to_string(r.position));
        value.init(v.gobj());
    }
    else if (column == columns.booksOwned.index() or column == columns.booksNew.index() or column == columns.booksWritten.index()
            or column == columns.lastBookAge.index()) {
        Glib::Value<size_t> v;
        v.init(v.value_type());
        v.set(  column == columns.booksOwned.index() ? r.library.size() :
                column == columns.booksNew.index() ? r.newBooks.size() :
                column == columns.booksWritten.index() ? r.wrote.size() :
                r.wrote.empty() ? state_.t : state_.books.at(r.wrote.back()).age
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
    if (sort_column_id == columns.id.index() || sort_column_id == DEFAULT_SORT_COLUMN_ID)
        compare = ascending ? less_id : greater_id;
    else if (sort_column_id == columns.posX.index())
        compare = ascending ? less_posX : greater_posX;
    else if (sort_column_id == columns.posY.index())
        compare = ascending ? less_posY : greater_posY;
    else if (sort_column_id == columns.posstr.index())
        compare = ascending ? less_posstr : greater_posstr;
    else if (sort_column_id == columns.u.index())
        compare = ascending ? less_uCurr : greater_uCurr;
    else if (sort_column_id == columns.uLifetime.index())
        compare = ascending ? less_uLife : greater_uLife;
    else if (sort_column_id == columns.booksOwned.index())
        compare = ascending ? less_booksOwned : greater_booksOwned;
    else if (sort_column_id == columns.booksNew.index())
        compare = ascending ? less_booksNew : greater_booksNew;
    else if (sort_column_id == columns.booksWritten.index())
        compare = ascending ? less_booksWritten : greater_booksWritten;
    else if (sort_column_id == columns.lastBookAge.index())
        // Can't handle like the above: less/greater_copies() is not static
        compare = std::bind(ascending ? &ReaderStore::less_lastBookAge : &ReaderStore::greater_lastBookAge, this, _1, _2);

    sort_members(compare, sort_column_id, order);
}

bool ReaderStore::less_id(const ReaderState &a, const ReaderState &b) {
    return a.id < b.id;
}
bool ReaderStore::greater_id(const ReaderState &a, const ReaderState &b) {
    return a.id > b.id;
}
bool ReaderStore::less_posX(const ReaderState &a, const ReaderState &b) {
    return a.position[0] < b.position[0];
}
bool ReaderStore::greater_posX(const ReaderState &a, const ReaderState &b) {
    return a.position[0] > b.position[0];
}
bool ReaderStore::less_posY(const ReaderState &a, const ReaderState &b) {
    return a.position[1] < b.position[1];
}
bool ReaderStore::greater_posY(const ReaderState &a, const ReaderState &b) {
    return a.position[1] > b.position[1];
}
// First x, then y for ties
bool ReaderStore::less_posstr(const ReaderState &a, const ReaderState &b) {
    auto ax = a.position[0], bx = b.position[0];
    return ax == bx ? a.position[1] < b.position[1] : ax < bx;
}
bool ReaderStore::greater_posstr(const ReaderState &a, const ReaderState &b) {
    auto ax = a.position[0], bx = b.position[0];
    return ax == bx ? a.position[1] > b.position[1] : ax > bx;
}
bool ReaderStore::less_uCurr(const ReaderState &a, const ReaderState &b) {
    return a.u < b.u;
}
bool ReaderStore::greater_uCurr(const ReaderState &a, const ReaderState &b) {
    return a.u > b.u;
}
bool ReaderStore::less_uLife(const ReaderState &a, const ReaderState &b) {
    return a.uLifetime < b.uLifetime;
}
bool ReaderStore::greater_uLife(const ReaderState &a, const ReaderState &b) {
    return a.uLifetime > b.uLifetime;
}
bool ReaderStore::less_booksOwned(const ReaderState &a, const ReaderState &b) {
    return a.library.size() < b.library.size();
}
bool ReaderStore::greater_booksOwned(const ReaderState &a, const ReaderState &b) {
    return a.library.size() > b.library.size();
}
bool ReaderStore::less_booksNew(const ReaderState &a, const ReaderState &b) {
    return a.newBooks.size() < b.newBooks.size();
}
bool ReaderStore::greater_booksNew(const ReaderState &a, const ReaderState &b) {
    return a.newBooks.size() > b.newBooks.size();
}
bool ReaderStore::less_booksWritten(const ReaderState &a, const ReaderState &b) {
    return a.wrote.size() < b.wrote.size();
}
bool ReaderStore::greater_booksWritten(const ReaderState &a, const ReaderState &b) {
    return a.wrote.size() > b.wrote.size();
}
bool ReaderStore::less_lastBookAge(const ReaderState &a, const ReaderState &b) const {
    return (a.wrote.empty() ? std::numeric_limits<unsigned long>::max() : state_.books.at(a.wrote.back()).age)
         < (b.wrote.empty() ? std::numeric_limits<unsigned long>::max() : state_.books.at(b.wrote.back()).age);
}
bool ReaderStore::greater_lastBookAge(const ReaderState &a, const ReaderState &b) const {
    return (a.wrote.empty() ? std::numeric_limits<unsigned long>::max() : state_.books.at(a.wrote.back()).age)
         > (b.wrote.empty() ? std::numeric_limits<unsigned long>::max() : state_.books.at(b.wrote.back()).age);
}

void ReaderStore::appendColumnsTo(Gtk::TreeView &v) const {
    appendCol(v, "ID", columns.id, 100);
    appendCol(v, "Position", columns.posstr, 150);
    appendCol(v, "Utility", columns.u, 100);
    appendCol(v, "Life Util.", columns.uLifetime, 100);
    appendCol(v, "Books", columns.booksOwned, 100);
    appendCol(v, "# New", columns.booksNew, 100);
    appendCol(v, "# Written", columns.booksWritten, 100);
    appendCol(v, "Last wrote", columns.lastBookAge, 100);
}

}}
