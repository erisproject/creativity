#pragma once
#include <gtkmm/treemodel.h>
#include <gtkmm/treeview.h>
#include <gtkmm/liststore.h>
#include <eris/SharedMember.hpp>

namespace creativity {

// Predeclarations:
class Book;

namespace gui {

/** Gtk::TreeModel::ColumnRecord subclass for handling Book information in the list of books in
 * the main GUI window.
 *
 * \sa GUI::AuthoredBookCols for a subclass of this class used for handling all Book information in
 * the reader/author dialog window, which works the same way but hides the author column.
 */
class BookCols : public Gtk::TreeModel::ColumnRecord {
    public:
        /** Creates an unpopulated column record */
        BookCols();
        /** Takes a Gtk::TreeView and adds this object's columns to it. */
        virtual void appendColumnsTo(Gtk::TreeView &v) const;
        /** Takes a list store and a reader to add a row using the reader's information to the list
         * store.
         */
        virtual void appendRow(Glib::RefPtr<Gtk::ListStore> &ls, const eris::SharedMember<Book> &book) const;

        Gtk::TreeModelColumn<unsigned long> id, author, age, sales, copies, lifetime;
        Gtk::TreeModelColumn<double> quality, price, revenue;
        Gtk::TreeModelColumn<bool> market;
        Gtk::TreeModelColumn<std::string> position;
};

} }
