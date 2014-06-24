#pragma once
#include <gtkmm/treemodel.h>
#include <gtkmm/treeview.h>
#include <gtkmm/liststore.h>
#include <eris/SharedMember.hpp>

namespace creativity {

// Predeclarations:
class Reader;

namespace gui {

/** Gtk::TreeModel::ColumnRecord subclass for handling Reader information in the list of readers in
 * the main GUI window.
 */
class ReaderCols : public Gtk::TreeModel::ColumnRecord {
    public:
        /** Creates an unpopulated column record */
        ReaderCols();
        /** Takes a Gtk::TreeView and adds this object's columns to it. */
        void appendColumnsTo(Gtk::TreeView &v) const;
        /** Takes a list store and a reader to add a row using the reader's information to the list
         * store.
         */
        void appendRow(Glib::RefPtr<Gtk::ListStore> &ls, const eris::SharedMember<Reader> &reader) const;

        Gtk::TreeModelColumn<unsigned long> id, ///< column for reader id
            books_owned, ///< column containing number of books owned
            books_new, ///< column containing number of books obtained in the previous period
            books_written, ///< column containing number of books authored by the reader
            book_latest_age; ///< column containing the age of the most recently authored book
        Gtk::TreeModelColumn<std::string> position; ///< column containing the reader position
        Gtk::TreeModelColumn<double> utility, ///< column containing the reader's previous period utility
            u_life; ///< column containing the reader's cumulative lifetime utility
};

} }
