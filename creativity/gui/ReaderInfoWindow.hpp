#pragma once
#include "creativity/gui/InfoWindow.hpp"

namespace creativity { namespace gui {

/** InfoWindow subclass for displaying the details of a reader/author.
 */
class ReaderInfoWindow : public InfoWindow {
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
        ReaderInfoWindow(std::shared_ptr<const state::State> state, std::shared_ptr<Gtk::Window> main_window,
                eris::eris_id_t reader_id, std::function<void(eris::eris_id_t)> open_info_dialog);

        /** Refresh the information in the dialog using the given simulation state. */
        virtual void refresh(std::shared_ptr<const state::State> state) override;

    private:
        std::function<void(eris::eris_id_t)> open_info_dialog_;

        Glib::RefPtr<BookStore> bk_authored_model_;
        Glib::RefPtr<LibraryStore> bk_library_model_;
        Gtk::TreeView bk_authored_tree_, bk_library_tree_;
};

}}
