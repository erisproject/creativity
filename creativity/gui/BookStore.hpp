#pragma once
#include "creativity/gui/MemberStore.hpp"

namespace creativity { namespace gui {

/** Gtk::TreeModel::ColumnRecord subclass for handling Book information in the list of books in
 * the main GUI window, and the list of authored books on the reader dialog.
 */
class BookStore : public MemberStore<state::BookState>, Glib::Object {
    public:
        BookStore() = delete;

        /** Interface class between a simulation's Books and a Gtk::TreeView.  This internally
         * stores a vector of Books which can be updated (if needed) by calling the update method.
         *
         * This exposes 9 columns:
         * - ID
         * - x position
         * - y position
         * - position string (made from x position and y position)
         * - author ID
         * - quality
         * - price (NaN if not on market)
         * - revenue
         * - revenueLifetime
         * - market (true if on market, false otherwise)
         * - age
         * - sales
         * - salesLifetime
         * - copies
         * - lifetime (# periods on market)
         *
         * \param state the simulation state object containing the books to list
         * \param author the author whose books to list.  Omit to list all state books.
         */
        static Glib::RefPtr<BookStore> create(const state::State &state, eris::eris_id_t author = 0);

        /** ColumnRecord object for a BookStore.  This object contains the columns for this Book
         * model.  This should not be used directly, but rather accessed via the public `columns`
         * member.
         */
        class ColRec : public Gtk::TreeModel::ColumnRecord {
            public:
                Gtk::TreeModelColumn<eris::eris_id_t> id; ///< Book ID
                Gtk::TreeModelColumn<eris::eris_id_t> author; ///< Author ID
                Gtk::TreeModelColumn<double> posX; ///< x coordinate of the book
                Gtk::TreeModelColumn<double> posY; ///< y coordinate of the book
                Gtk::TreeModelColumn<std::string> posstr; ///< position of the book as a string such as `(-7.16,0.440)`
                Gtk::TreeModelColumn<double> quality; ///< quality parameter of the book (the mean of realized quality draws)
                Gtk::TreeModelColumn<bool> market; ///< True if the book is currently on the market
                Gtk::TreeModelColumn<double> price; ///< price of the book, or NaN if the book is not on the market
                Gtk::TreeModelColumn<double> revenue; ///< revenue of the book in the current period
                Gtk::TreeModelColumn<double> revenueLifetime; ///< Cumulative revenue of the book since its creation
                Gtk::TreeModelColumn<size_t> age; ///< Age of the book in simulation periods since it was written
                Gtk::TreeModelColumn<size_t> sales; ///< Copies sold in the current period
                Gtk::TreeModelColumn<size_t> salesLifetime; ///< Lifetime copies sold
                /** Copies that exist in the simulation.  This is at least one larger than the
                 * number of lifetime sales because the author has a copy (which wasn't a sale); if
                 * there is non-sale piracy, this value could be much greater than lifetime sales.
                 */
                Gtk::TreeModelColumn<size_t> copies;
                Gtk::TreeModelColumn<size_t> lifetime; ///< Number of periods the book has been or was on the market

            private:
                ColRec() {
                    add(id); add(author); add(market); add(posX); add(posY); add(posstr);
                    add(quality); add(price); add(revenue); add(revenueLifetime);
                    add(age); add(sales); add(salesLifetime); add(copies); add(lifetime);
                }
                friend class BookStore;
        };

        /** The columns of this BookStore.  For example, to access the price column, use
         * `bookstore.columns.price`.
         */
        ColRec columns;

        /** Takes a Gtk::TreeView and adds this object's columns to it. */
        void appendColumnsTo(Gtk::TreeView &v) const;

    protected:
        /// Protected constructor; this object should be constructed using create().
        BookStore(const state::State &state, eris::eris_id_t author = 0);

        /** Returns `obj.columns.size()`, the number of book model columns. */
        virtual int get_n_columns_vfunc() const override;

        /** Returns the column type of the given position.  See the list of virtual columns in the
         * class description.
         *
         * \sa BookStore::ColRec
         */
        virtual GType get_column_type_vfunc(int index) const override;

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
        // If non-zero, this is an author-specific list; otherwise it's a global list
        const eris::eris_id_t author_ = 0;

        // The various comparison functions; one of these gets passed to std::stable_sort.
#define LESS_GREATER_METHODS(col) \
        static bool less_##col(const state::BookState &a, const state::BookState &b); \
        static bool greater_##col(const state::BookState &a, const state::BookState &b);
        LESS_GREATER_METHODS(id)
        LESS_GREATER_METHODS(author)
        LESS_GREATER_METHODS(market)
        LESS_GREATER_METHODS(posX)
        LESS_GREATER_METHODS(posY)
        LESS_GREATER_METHODS(quality)
        LESS_GREATER_METHODS(price)
        LESS_GREATER_METHODS(revenue)
        LESS_GREATER_METHODS(revenueLifetime)
        LESS_GREATER_METHODS(posstr)
        LESS_GREATER_METHODS(age)
        LESS_GREATER_METHODS(sales)
        LESS_GREATER_METHODS(salesLifetime)
        LESS_GREATER_METHODS(lifetime)
        LESS_GREATER_METHODS(copies)
#undef LESS_GREATER_METHODS
};

}}
