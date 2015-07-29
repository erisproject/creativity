#pragma once
#include "creativity/gui/InfoWindow.hpp"

namespace creativity { namespace gui {

/** InfoWindow subclass for displaying the details of a book.
 */
class BookInfoWindow : public InfoWindow {
    public:
        /** Constructs a new InfoWindow displaying book information.
         *
         * \param state the state at which to display the information for the initial refresh()
         * \param main_window the main window to which this dialog will be attached
         * \param book_id the id of the book to display.  If the book does not exist in the
         * passed-in state, its values will be set to "N/A".  (This is allowed because the book
         * window can be refresh()ed to another state where the book may or may not exist).
         */
        BookInfoWindow(std::shared_ptr<const state::State> state, std::shared_ptr<Gtk::Window> main_window,
                eris::eris_id_t book_id);

        /** Refresh the information in the dialog using the given simulation state. */
        virtual void refresh(std::shared_ptr<const state::State> state) override;
        
};

}}
