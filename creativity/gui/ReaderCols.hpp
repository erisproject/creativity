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

        Gtk::TreeModelColumn<unsigned long> id, books_owned, books_new, books_written, book_latest_age;
        Gtk::TreeModelColumn<std::string> position;
        Gtk::TreeModelColumn<double> utility, u_life;
};

} }
