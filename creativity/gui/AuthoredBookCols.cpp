#include "creativity/gui/AuthoredBookCols.hpp"
#include "creativity/Reader.hpp"

using namespace eris;

namespace creativity { namespace gui {

void AuthoredBookCols::appendColumnsTo(Gtk::TreeView &v) const {
    BookCols::appendColumnsTo(v);
    // NB: if BookCols column position of "Author" changes, change this too:
    v.remove_column(*v.get_column(2));
}

} }
