#pragma once
#include "creativity/gui/BookStore.hpp"
#include "creativity/state/BookState.hpp"
#include <eris/types.hpp>
#include <glibmm/object.h>
#include <glibmm/refptr.h>
#include <gtkmm/enums.h>
#include <gtkmm/treemodel.h>
#include <gtkmm/treemodelcolumn.h>
#include <memory>

namespace Glib { class ValueBase; }
namespace Gtk { class TreeView; }
namespace creativity { namespace state { class State; } }

namespace creativity { namespace gui {

/** Class extending BookStore that lists a reader's library books (excluding self-written works); it
 * includes all the base book columns plus:
 * - perceived quality (double)
 * - pirated (boolean)
 * - bought (boolean)
 * - acquired (simulation period)
 */
class LibraryStore : public BookStore, public virtual Glib::Object {
    public:
        /** Creates a LibraryStore (cast to a BookStore) for the given reader. */
        static Glib::RefPtr<LibraryStore> create(std::shared_ptr<const state::State> state, eris::eris_id_t reader);

        /** BookStore::ColRec extension to add reader-specific fields */
        class ColRec : public BookStore::ColRec {
            public:
                Gtk::TreeModelColumn<double> reader_quality; ///< This reader's quality draw for the book

                Gtk::TreeModelColumn<bool>
                    reader_pirated, ///< True if the book was obtained by piracy
                    reader_purchased; ///< True if the book was purchased

                Gtk::TreeModelColumn<unsigned int>
                    reader_acquired; ///< Simulation period when this copy of the book was acquired

            protected:
                ColRec() : BookStore::ColRec() {
                    add(reader_quality);
                    add(reader_pirated);
                    add(reader_purchased);
                    add(reader_acquired);
                }
                friend class LibraryStore;
        };

        /** Takes a Gtk::TreeView and adds this object's columns to it. */
        virtual void appendColumnsTo(Gtk::TreeView &v) const override;

    protected:
        /// Protected constructor; this object should be constructed using create().
        LibraryStore(std::shared_ptr<const state::State> &&state, eris::eris_id_t reader);

        /** Accesses a column value.
         *
         * \param iter a valid iterator referencing the row to access
         * \param column the index of the column to access
         * \param value a Glib::Value<TYPE> object (where TYPE is the appropriate type for the
         * requested `column`) in which to store the value.
         */
        virtual void get_value_vfunc(const iterator &iter, int column, Glib::ValueBase &value) const override;

        // TreeSortable overrides:

        /** Sets the model sort column and sort order.  If the sort_column and order differ from the
         * current values, the model data is resorted.  If resorting occurs, the sort_column_changed
         * and rows_reordered signals will fire.
         *
         * This uses a stable sort: elements that are equal under the new sort column will preserve
         * their current ordering.  Note that this means that sorting by the current column but in
         * the opposite order will *not* reverse the ordering of equal-value elements.
         *
         * \param sort_column_id the index of the new sort column
         * \param order the new sort order (Gtk::SORT_ASCENDING or Gtk::SORT_DESCENDING).
         */
        virtual void set_sort_column_id_vfunc(int sort_column_id, Gtk::SortType order) override;

    private:
        // The reader whose library this Store represents
        const eris::eris_id_t reader_;

        // The various comparison functions; one of these gets passed to std::stable_sort.
#define LESS_GREATER_METHODS(col) \
        bool less_##col(const state::BookState &a, const state::BookState &b); \
        bool greater_##col(const state::BookState &a, const state::BookState &b);
        LESS_GREATER_METHODS(reader_quality)
        LESS_GREATER_METHODS(reader_purchased)
        LESS_GREATER_METHODS(reader_pirated)
        LESS_GREATER_METHODS(reader_acquired)
#undef LESS_GREATER_METHODS
};

}}
