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

using namespace Eigen;

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
        InfoWindow(std::shared_ptr<const state::State> state, std::shared_ptr<Gtk::Window> main_window,
                eris::eris_id_t reader_id, std::function<void(eris::eris_id_t)> open_info_dialog);

        /** Constructs a new InfoWindow displaying book information.
         *
         * \param state the state at which to display the information for the initial refresh()
         * \param main_window the main window to which this dialog will be attached
         * \param book_id the id of the book to display.  If the book does not exist in the
         * passed-in state, its values will be set to "N/A".  (This is allowed because the book
         * window can be refresh()ed to another state where the book may or may not exist).
         */
        InfoWindow(std::shared_ptr<const state::State> state, std::shared_ptr<Gtk::Window> main_window,
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
        void refresh(std::shared_ptr<const state::State> state);

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
        /** Updates a matrix set up with matrix_at with the values of the given matrix.  If the
         * optional lower_triangle parameter is given as true, only the lower diagonal of the matrix
         * is displayed; the upper-triangle values are set to blanks.
         */
        void updateMatrix(const std::string &code, const Ref<const MatrixXd> &m, bool lower_triangle = false);

    private:
        std::unordered_map<std::string, std::pair<Gtk::Label, Gtk::Label>> fields_;
        std::unordered_map<std::string, std::vector<std::unique_ptr<Gtk::Label>>> matrix_;
        std::list<Gtk::Grid> grids_;
        std::list<Gtk::Label> labels_;
        std::list<Gtk::Notebook> nbs_;
        std::list<Gtk::ScrolledWindow> swins_;

        // The global label position, used during setup.
        int gpos_ = 0;

        // Initialization helper functions:

        // Create a new grid, store it, return a reference to it (but don't add to any element)
        Gtk::Grid& new_grid();
        // Create a new grid and adds it as the content of new tab "title" on `notebook`, returns a
        // reference to it.
        Gtk::Grid& new_tab_grid(Gtk::Notebook &notebook, const std::string &title);
        // Creates a "value: data" label pair with labels stored in fields_ that can be updated
        // later.  The pair of labels is located at (`row`,`col`) and (`row`,`col+1`).
        void data_at(Gtk::Grid &grid, const std::string &code, const std::string &value_name, int row, int col);
        // Calls data_at to put the label at (`gpos_`, 0), then increments gpos_.
        void data_append(Gtk::Grid &grid, const std::string &code, const std::string &value_name);
        // Creates a label at the given location, taking up the given number of rows and columns,
        // with horizontal alignment as given (0-1).  Intended for static values.
        void label_at(Gtk::Grid &grid, const std::string &label, double alignment, int row, int col, int width, int height);
        // Creates a pair of single-cell labels at the given (row,col) and (row,col+1).  The first
        // is right-aligned, the second is left-aligned.
        void labels_at(Gtk::Grid &grid, const std::string &label1, const std::string &label2, int row, int col);
        // Creates a pair of static labels at the current (gpos_) row, first two columns.
        void labels_append(Gtk::Grid &grid, const std::string &label1, const std::string &label2);
        // Creates a static, wrapping, multi-column label, at the current (gpos_) row.
        void comment_append(Gtk::Grid &grid, const std::string &comment, int cols = 2);
        // Creates labels for a matrix of the given size in the grid with [0,0] element at
        // (`row+1`,`col+1`), storing the labels in row-major order in matrix_.  updateMatrix() can
        // later be called to update the label values to a given matrix.
        // (`row`,`col`) gets the given label, and the remainder of the row and column get position
        // labels.
        void matrix_at(Gtk::Grid &grid, const std::string &code, const std::string &name, int row, int col, int nrows, int ncols);

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
