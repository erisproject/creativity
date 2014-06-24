#include "creativity/gui/ReaderCols.hpp"
#include "creativity/gui/GUI.hpp"
#include "creativity/Reader.hpp"

using namespace eris;

namespace creativity { namespace gui {

ReaderCols::ReaderCols() {
    add(id);
    add(position);
    add(utility);
    add(u_life);
    add(books_owned);
    add(books_new);
    add(books_written);
    add(book_latest_age);
}

void ReaderCols::appendColumnsTo(Gtk::TreeView &v) const {
#define ADDCOL(label, col, sortable) \
    v.append_column(label, col);\
    if (sortable) v.get_column(v.get_n_columns()-1)->set_sort_column(col);
    ADDCOL("ID", id, 1);
    ADDCOL("Position", position, 0);
    ADDCOL("Utility", utility, 1);
    ADDCOL("Life Util.", u_life, 1);
    ADDCOL("Books", books_owned, 1);
    ADDCOL("# New", books_new, 1);
    ADDCOL("# Written", books_written, 1);
    ADDCOL("Last wrote", book_latest_age, 1);
}

void ReaderCols::appendRow(Glib::RefPtr<Gtk::ListStore> &ls, const SharedMember<Reader> &reader) const {
    auto row = *ls->append();
    row[id] = reader->id();
    row[position] = GUI::pos_to_string(reader->position());
    row[utility] = reader->u();
    row[u_life] = reader->uLifetime();
    row[books_owned] = reader->library().size();
    row[books_new] = reader->newBooks().size();
    row[books_written] = reader->wrote().size();
    if (not reader->wrote().empty())
        row[book_latest_age] = reader->wrote().back()->age();
}

} }
