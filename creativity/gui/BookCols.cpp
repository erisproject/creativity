#include "creativity/gui/BookCols.hpp"
#include "creativity/gui/GUI.hpp"
#include "creativity/Book.hpp"
#include "creativity/Reader.hpp"

using namespace eris;

namespace creativity { namespace gui {

BookCols::BookCols() {
    add(id);
    add(position);
    add(author);
    add(age);
    add(quality);
    add(market);
    add(price);
    add(revenue);
    add(sales);
    add(copies);
    add(lifetime);
}

void BookCols::appendColumnsTo(Gtk::TreeView &v) const {
#define ADDCOL(label, col, sortable) \
    v.append_column(label, col);\
    if (sortable) v.get_column(v.get_n_columns()-1)->set_sort_column(col);
    // NB: If changing order, adding or removing, ensure that either author's numeric position
    // doesn't change, or else that AuthoredBookCols gets updated with the new position.
    ADDCOL("ID", id, 1);
    ADDCOL("Position", position, 0);
    ADDCOL("Author", author, 1); // Position important!
    ADDCOL("Age", age, 1);
    ADDCOL("Quality", quality, 1);
    ADDCOL("Mkt?", market, 1);
    ADDCOL("Price", price, 1);
    ADDCOL("Rev.", revenue, 1);
    ADDCOL("Sales", sales, 1);
    ADDCOL("Copies", copies, 1);
    ADDCOL("Life", lifetime, 1);
}

void BookCols::appendRow(Glib::RefPtr<Gtk::ListStore> &ls, const SharedMember<Book> &book) const {
    auto row = *ls->append();
    row[id] = book->id();
    row[position] = GUI::pos_to_string(book->position());
    row[author] = book->author()->id();
    row[age] = book->age();
    row[quality] = book->quality();
    row[market] = book->hasMarket();
    row[price] = book->hasMarket() ? book->market()->price() : std::numeric_limits<double>::quiet_NaN();
    row[revenue] = book->lifeRevenue();
    row[sales] = book->lifeSales();
    row[copies] = book->copies();
    row[lifetime] = book->marketPeriods();
}

} }
