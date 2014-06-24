#pragma once
#include "creativity/gui/BookCols.hpp"

namespace creativity { namespace gui {

/** Subclass of BookCols that hides the author column, for use in an author's details dialog window.
 *
 * \sa creativity::GUI::BookCols
 */
class AuthoredBookCols : public BookCols {
    public:
        /** Enhances BookCols::appendColumnsTo by removing the Author column. */
        void appendColumnsTo(Gtk::TreeView &v) const override;
};

} }
