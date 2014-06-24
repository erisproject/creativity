#pragma once
#include <eris/Simulation.hpp>
#include <eris/noncopyable.hpp>
#include <eris/Position.hpp>
#include <gtkmm/window.h>
#include <gtkmm/notebook.h>
#include <gtkmm/grid.h>
#include <gtkmm/label.h>
#include <gtkmm/button.h>
#include <gtkmm/liststore.h>
#include <gtkmm/treeview.h>
#include <gtkmm/scrolledwindow.h>
#include <mutex>
#include <memory>
#include "creativity/gui/AuthoredBookCols.hpp"

namespace creativity {

// forward declarations
class Reader;
class Book;

namespace gui {

class GUI;

/** Gtk dialog for showing reader or book info. */
class InfoWindow : public Gtk::Window {
    public:
        /** Constructs a new InfoWindow of the given size displaying reader information.
         */
        InfoWindow(eris::SharedMember<Reader> reader, std::shared_ptr<Gtk::Window> main_window);

        /** Constructs a new InfoWindow of the given size displaying book information.
         */
        InfoWindow(eris::SharedMember<Book> book, std::shared_ptr<Gtk::Window> main_window);

        /** Refresh the information in the dialog. */
        void refresh();

        /** The reader about which this dialog is displaying details, or a null SharedMember is this
         * is a book details dialog.
         */
        eris::SharedMember<Reader> reader;

        /** The book about which this dialog is displaying details, or a null SharedMember is this
         * is a reader details dialog.
         */
        eris::SharedMember<Book> book;
    protected:
        /// Updates a single value Gtk::Label text with the given string
        void updateValue(const std::string &code, const std::string &val);
        /// Updates a single value Gtk::Label text with the given unsigned long
        void updateValue(const std::string &code, unsigned long val);
        /// Updates a single value Gtk::Label text with the given long
        void updateValue(const std::string &code, long val);
        /// Updates a single value Gtk::Label text with the given double
        void updateValue(const std::string &code, double val);

    private:
        std::unordered_map<std::string, std::pair<Gtk::Label, Gtk::Label>> fields_;
        std::list<Gtk::Grid> grids_;
        std::list<Gtk::Label> labels_;
        std::list<Gtk::Notebook> nbs_;
        std::list<Gtk::ScrolledWindow> swins_;

        AuthoredBookCols abc_;
        Glib::RefPtr<Gtk::ListStore> bk_list_;
        Gtk::TreeView bk_tree_;
};

} }
