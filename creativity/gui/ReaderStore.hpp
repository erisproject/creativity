#pragma once
#include <eris/types.hpp>
#include "creativity/gui/MemberStore.hpp"
#include "creativity/state/State.hpp"
#include <memory>
#include <cstddef>
#include <string>
#include <glibmm/object.h>
#include <glibmm/refptr.h>
#include <gtkmm/enums.h>
#include <gtkmm/treemodel.h>
#include <gtkmm/treemodelcolumn.h>
#include <gtkmm/treepath.h>
namespace Glib { class ValueBase; }
namespace Gtk { class TreeView; }

namespace creativity { namespace gui {

/** Gtk::TreeModel::ColumnRecord subclass for handling Reader information in the list of readers in
 * the main GUI window.
 */
class ReaderStore : public MemberStore<state::ReaderState>, private virtual Glib::Object {
    public:
        /** Interface class between a simulation state's Reader information and a Gtk::TreeView.
         *
         * This exposes various reader-specific columns; see ColRec for details.
         */
        static Glib::RefPtr<ReaderStore> create(std::shared_ptr<const state::State> state);

        /** ColumnRecord object for a ReaderStore.  This object contains the columns for this Book
         * model.  This should not be used directly, but rather accessed via the public `columns`
         * member.
         */
        class ColRec : public Gtk::TreeModel::ColumnRecord {
            public:
                Gtk::TreeModelColumn<eris::eris_id_t> id; ///< ID of the reader
                Gtk::TreeModelColumn<double> pos_x, ///< X coordinate of the reader
                    pos_y, ///< y coordinate of the reader
                    u, ///< Current period utility of the reader
                    u_lifetime, ///< Cumulative lifetime utility of the reader
                    creation_shape, ///< The reader's creation_shape parameter
                    creation_scale; ///< The reader's creation_scale parameter
                Gtk::TreeModelColumn<std::string> pos_str; ///< position of the book as a string such as `(-7.16,0.440)`
                Gtk::TreeModelColumn<size_t> books_owned, ///< Number of books possessed by this reader
                    books_purchased, ///< Number of books purchased by this reader
                    books_pirated, ///< Number of books obtained via piracy by this reader
                    books_new, ///< Number of books acquired in the current period
                    books_new_purchased, ///< Number of books acquired by purchasing in the current period
                    books_new_pirated, ///< Number of books acquired by piracy in the current period
                    books_written, ///< Number of books authored by this reader
                    last_book_age, ///< Age of the most recently written book, or simulation age if `booksWritten == 0`
                    num_friends; ///< Number of friends of this reader

            private:
                ColRec() {
                    add(id); add(pos_x); add(pos_y); add(pos_str); add(u); add(u_lifetime);
                    add(creation_shape); add(creation_scale); add(num_friends);
                    add(books_owned); add(books_purchased); add(books_pirated);
                    add(books_new); add(books_new_purchased); add(books_new_pirated);
                    add(books_written); add(last_book_age);
                }
                friend class ReaderStore;
        };

        /** The columns of this ReaderStore.  For example, to access the utility column, use
         * `readerstore.columns.u`.
         */
        ColRec columns;

        /** Takes a Gtk::TreeView and adds this object's columns to it. */
        void appendColumnsTo(Gtk::TreeView &v) const override;

    protected:
        /// Protected constructor; this object should be constructed using create().
        ReaderStore(std::shared_ptr<const state::State> &&state);

        /** Returns the column type of the given position.  This is typically invoked via
         * get_column_type, itself given a column member of the `.columns` ColRec object.
         *
         * \sa ReaderStore::ColRec
         */
        virtual GType get_column_type_vfunc(int index) const override;

        /** Returns `obj.columns.size()`, the number of model columns. */
        virtual int get_n_columns_vfunc() const override;

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
        // The various comparison functions; one of these gets passed to std::stable_sort.
#define LESS_GREATER_METHODS(col) \
        static bool less_##col(const state::ReaderState &a, const state::ReaderState &b); \
        static bool greater_##col(const state::ReaderState &a, const state::ReaderState &b);
        LESS_GREATER_METHODS(id)
        LESS_GREATER_METHODS(pos_x)
        LESS_GREATER_METHODS(pos_y)
        LESS_GREATER_METHODS(pos_str)
        LESS_GREATER_METHODS(u)
        LESS_GREATER_METHODS(u_lifetime)
        LESS_GREATER_METHODS(creation_shape)
        LESS_GREATER_METHODS(creation_scale)
        LESS_GREATER_METHODS(books_owned)
        LESS_GREATER_METHODS(books_purchased)
        LESS_GREATER_METHODS(books_pirated)
        LESS_GREATER_METHODS(books_new)
        LESS_GREATER_METHODS(books_new_purchased)
        LESS_GREATER_METHODS(books_new_pirated)
        LESS_GREATER_METHODS(books_written)
        LESS_GREATER_METHODS(num_friends)
#undef LESS_GREATER_METHODS
        bool less_last_book_age(const state::ReaderState &a, const state::ReaderState &b) const;
        bool greater_last_book_age(const state::ReaderState &a, const state::ReaderState &b) const;

};

}}
