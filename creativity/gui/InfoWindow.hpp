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
#include <Eigen/Core>
#include "creativity/gui/BookStore.hpp"
#include "creativity/gui/LibraryStore.hpp"
#include "creativity/state/State.hpp"

namespace creativity {

// forward declarations
class Reader;
class Book;

namespace gui {

class GUI;

/** Gtk dialog for showing reader or book info.  Base class; see ReaderInfoWindow and
 * BookInfoWindow. */
class InfoWindow : public Gtk::Window {
    public:
        /// No default constructor
        InfoWindow() = delete;

        /** The eris_id of the simulation object about which this dialog is displaying details.
         */
        const eris::eris_id_t id;

        /** Refresh the information in the dialog using the given simulation state. */
        virtual void refresh(std::shared_ptr<const state::State> state) = 0;

    protected:
        /** Creates a new InfoWindow, setting up the dialog window.
         *
         * \param main_window the parent Gtk::Window this window attaches to
         * \param id the eris id of the displayed element
         * \param title the window title
         */
        InfoWindow(std::shared_ptr<Gtk::Window> main_window, eris::eris_id_t id, const std::string &title);

        /** Updates a single value Gtk::Label text with the given value by passing it to
         * std::to_string.
         */
        template <class T, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
        void updateValue(const std::string &code, const T &val) {
            updateValue(code, std::to_string(val));
        }

        /// Updates a single value Gtk::Label text with the given string
        void updateValue(const std::string &code, const std::string &val);

        /** Updates a matrix set up with matrix_at with the values of the given matrix.  If the
         * optional lower_triangle parameter is given as true, only the lower diagonal of the matrix
         * is displayed; the upper-triangle values are set to blanks.
         */
        void updateMatrix(const std::string &code, const Eigen::Ref<const Eigen::MatrixXd> &m, bool lower_triangle = false);

        /** Clears all values of a matrix set up with matrix_at, setting all values to the given
         * value (defaulting to an empty string).
         */
        void clearMatrix(const std::string &code, const std::string &value = "");

        /// Fields containing a pair of labels: description and value
        std::unordered_map<std::string, std::pair<Gtk::Label, Gtk::Label>> fields_;
        /// Field containing a description and a matrix of values.
        std::unordered_map<std::string, std::vector<std::unique_ptr<Gtk::Label>>> matrix_;
        /// Grid for the display
        std::list<Gtk::Grid> grids_;
        /// Stored labels that aren't fields (i.e. single labels)
        std::list<Gtk::Label> labels_;
        /// Notebooks for the tabs
        std::list<Gtk::Notebook> nbs_;
        /// Any scrolled windows used
        std::list<Gtk::ScrolledWindow> swins_;

        /// The global label position, used during setup.
        int gpos_ = 0;

        // Initialization helper functions:

        /// Create a new grid, store it, return a reference to it (but don't add to any element)
        Gtk::Grid& new_grid();
        /** Create a new grid and adds it as the content of new tab "title" on `notebook`, returns a
         * reference to it.
         */
        Gtk::Grid& new_tab_grid(Gtk::Notebook &notebook, const std::string &title);
        /* Creates a "value: data" label pair with labels stored in fields_ that can be updated
         * later.  The pair of labels is located at (`row`,`col`) and (`row`,`col+1`).
         */
        void data_at(Gtk::Grid &grid, const std::string &code, const std::string &value_name, int row, int col,
                double label_align = 1.0, double value_align = 0.0);
        /// Calls data_at to put the label at (`gpos_`, 0), then increments gpos_.
        void data_append(Gtk::Grid &grid, const std::string &code, const std::string &value_name,
                double label_align = 1.0, double value_align = 0.0);
        /** Creates a label at the given location, taking up the given number of rows and columns,
         * with horizontal alignment as given (0-1).  Intended for static values.
         */
        void label_at(Gtk::Grid &grid, const std::string &label, double alignment, int row, int col, int width, int height);
        /** Creates a pair of single-cell labels at the given (row,col) and (row,col+1).  The first
         * is right-aligned, the second is left-aligned.
         */
        void labels_at(Gtk::Grid &grid, const std::string &label1, const std::string &label2, int row, int col);
        /// Creates a pair of static labels at the current (gpos_) row, first two columns.
        void labels_append(Gtk::Grid &grid, const std::string &label1, const std::string &label2);
        /// Creates a static, wrapping, multi-column label, at the current (gpos_) row.
        void comment_append(Gtk::Grid &grid, const std::string &comment, int cols = 2);
        /** Creates labels for a matrix of the given size in the grid with [0,0] element at
         * (`row+1`,`col+1`), storing the labels in row-major order in matrix_.  updateMatrix() can
         * later be called to update the label values to a given matrix.  (`row`,`col`) gets the
         * given label, and the remainder of the row and column get position labels.
         */
        void matrix_at(Gtk::Grid &grid, const std::string &code, const std::string &name, int row, int col, int nrows, int ncols);

        // The time currently being shown (don't need to do anything in refresh() if this doesn't
        // change)
        unsigned long t_ = (unsigned long) -1;
        // On the initial refresh, don't replace the book model
        bool initial_refresh_ = true;
};

} }
