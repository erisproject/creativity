#pragma once
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
#include "creativity/gui/BookStore.hpp"
#include "creativity/state/State.hpp"

namespace creativity {

// forward declarations
class Reader;
class Book;

namespace gui {

class GUI;

/** Gtk dialog for showing reader or book info. */
class InfoWindow : public Gtk::Window {
    public:
        /** Constructs a new InfoWindow displaying reader information.
         *
         * \param state the state at which to display the information for the initial refresh()
         * \param main_window the main window to which this dialog will be attached
         * \param reader_id the reader id to display, or 0 to display a book
         * \param open_info_dialog a callable object that takes an ID and opens an info dialog (used
         * to open authored book dialogs).
         *
         * \throws std::out_of_range if the reader or book doesn't exist in the passed-in state.
         */
        InfoWindow(const state::State &state, std::shared_ptr<Gtk::Window> main_window,
                eris::eris_id_t reader_id, std::function<void(eris::eris_id_t)> open_info_dialog);

        /** Constructs a new InfoWindow displaying book information.
         *
         * \param state the state at which to display the information for the initial refresh()
         * \param main_window the main window to which this dialog will be attached
         * \param book_id the id of the book to display.  If the book does not exist in the
         * passed-in state, its values will be set to "N/A".  (This is allowed because the book
         * window can be refresh()ed to another state where the book may or may not exist).
         */
        InfoWindow(const state::State &state, std::shared_ptr<Gtk::Window> main_window,
                eris::eris_id_t book_id);

        /** The reader id about which this dialog is displaying details, or 0 if this is a book
         * details dialog.
         */
        const eris::eris_id_t reader;

        /** The book about which this dialog is displaying details, or 0 if this is a reader details
         * dialog.
         */
        const eris::eris_id_t book;

        /** Refresh the information in the dialog using the given simulation state. */
        void refresh(const state::State &state);

    protected:
        /// Sets up the dialog window properties, associates it with the parent, etc.
        void initWindow(Gtk::Window &parent);
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

        std::function<void(eris::eris_id_t)> open_info_dialog_;

        Glib::RefPtr<BookStore> bk_model_;
        Gtk::TreeView bk_tree_;
        
        // The time currently being shown (don't need to do anything in refresh() if this doesn't
        // change)
        unsigned long t_ = (unsigned long) -1;
        // On the initial refresh, don't replace the book model
        bool initial_refresh_ = true;
};

} }
