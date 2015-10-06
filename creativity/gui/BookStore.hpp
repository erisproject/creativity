#pragma once
#include "creativity/gui/MemberStore.hpp"
#include "creativity/state/BookState.hpp"
#include <eris/types.hpp>
#include <glibmm/object.h>
#include <glibmm/refptr.h>
#include <gtkmm/enums.h>
#include <gtkmm/treemodel.h>
#include <gtkmm/treemodelcolumn.h>
#include <gtkmm/treepath.h>
#include <memory>
#include <string>

namespace Glib { class ValueBase; }
namespace Gtk { class TreeView; }
namespace creativity { namespace state { class State; } }


namespace creativity { namespace gui {

/** Gtk::TreeModel::ColumnRecord subclass for handling Book information in the list of books in
 * the main GUI window, and the list of authored books on the reader dialog.
 */
class BookStore : public MemberStore<state::BookState>, public virtual Glib::Object {
    public:
        BookStore() = delete;

        /** Interface class between a simulation's Books and a Gtk::TreeView.  This internally
         * stores a vector of Books which can be updated (if needed) by calling the update method.
         *
         * This exposes the following columns:
         * - ID
         * - author ID
         * - x position
         * - y position
         * - position string (made from x position and y position)
         * - quality (reader-perceived quality may diff)
         * - price (NaN if not on market)
         * - revenue
         * - revenue_lifetime
         * - market_private (true if on private market, false otherwise)
         * - market_public (true if on public market, false otherwise)
         * - market_any (true if on any market, false otherwise)
         * - market_str (string: "Priv.", "Pub.", or "No")
         * - age (in simulation periods)
         * - creation date (i.e. simulation period)
         * - current sales
         * - lifetime private sales
         * - lifetime public sales
         * - lifetime sales (private + public)
         * - current pirated copies
         * - lifetime pirated copies
         * - lifetime copies (= lifetime pirated copies + lifetime sales)
         * - lifetime (# periods on market)
         *
         * \param state the simulation state object containing the books to list
         * \param author the author whose books to list.  Omit or set to 0 to list all books.
         */
        static Glib::RefPtr<BookStore> create(std::shared_ptr<const state::State> state, eris::eris_id_t author = 0);

        /** ColumnRecord object for a BookStore.  This object contains the columns for this Book
         * model.  This should not be used directly, but rather accessed via the public `columns`
         * member.
         */
        class ColRec : public Gtk::TreeModel::ColumnRecord {
            public:
                Gtk::TreeModelColumn<eris::eris_id_t>
                    id, ///< Book ID
                    author; ///< Author ID

                Gtk::TreeModelColumn<double>
                    pos_x, ///< x coordinate of the book
                    pos_y, ///< y coordinate of the book
                    quality, ///< quality parameter of the book (the mean of realized quality draws)
                    price, ///< price of the book, or NaN if the book is not on the market
                    revenue, ///< revenue of the book in the current period
                    revenue_lifetime; ///< Cumulative revenue of the book since its creation

                Gtk::TreeModelColumn<std::string>
                    pos_str, ///< position of the book as a string such as `(-7.16,0.440)`
                    market_str; ///< market type string ("Priv.", "Pub.", or "No")

                Gtk::TreeModelColumn<bool>
                    market_private, ///< True if on the market and that market is private
                    market_public, ///< True if on the market and that market is public
                    market_any; ///< True if the book is currently on either market

                Gtk::TreeModelColumn<unsigned int>
                    age, ///< Age of the book in simulation periods since it was written
                    created, ///< Age of the book in simulation periods since it was written
                    sales, ///< Copies sold in the current period
                    sales_lifetime_private, ///< Lifetime copies sold on the private market
                    sales_lifetime_public, ///< Lifetime copies sold on the public market
                    sales_lifetime, ///< Lifetime copies sold (both private/public)
                    pirated, ///< Copies sold in the current period
                    pirated_lifetime, ///< Lifetime copies sold
                    copies, ///< Copies created (sold or pirated) in the current period
                    copies_lifetime, ///< Lifetime copies created (sold or pirated)
                    lifetime_private; ///< Number of periods the book has been or was on the private market

            protected:
                ColRec() {
                    add(id); add(author);
                    add(pos_x); add(pos_y); add(quality); add(price); add(revenue); add(revenue_lifetime);
                    add(pos_str);
                    add(market_private); add(market_public); add(market_any); add(market_str);
                    add(age); add(created); add(sales); add(sales_lifetime_private); add(sales_lifetime_public);
                    add(sales_lifetime); add(pirated); add(pirated_lifetime);
                    add(copies); add(copies_lifetime); add(lifetime_private);
                }
                friend class BookStore;
        };

        /** The columns of this BookStore.  For example, to access the price column, use
         * `bookstore.columns->price`.
         */
        std::unique_ptr<ColRec> columns{new ColRec};

        /** Takes a Gtk::TreeView and adds this object's columns to it. */
        virtual void appendColumnsTo(Gtk::TreeView &v) const override;

    protected:
        /// Protected constructor; this object should be constructed using create().
        BookStore(std::shared_ptr<const state::State> &&state, eris::eris_id_t author);

        /// Protected constructor that overrides the default set of columns (for use by subclasses).
        BookStore(std::shared_ptr<const state::State> &&state, eris::eris_id_t author, std::unique_ptr<ColRec> &&cols);

        /** Protected constructor that overrides the default set of columns but has no author
         * parameter at all: members_ will not be populated (but should be populated by the
         * subclass).
         */
        BookStore(std::shared_ptr<const state::State> &&state, std::unique_ptr<ColRec> &&cols);

        /** Called during construction to populated members_ with either all books (if the author
         * constructor parameter is 0) or an author's written books (if the author constructor
         * parameter is non-zero).  This is not called in the subclass constructor that has no
         * author parameter (explicit or implicit) at all.
         */
        void initializeBooks();

        /** Returns `obj.columns->size()`, the number of book model columns. */
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
        // Filter by author, if non-zero; otherwise, this is a global book list (or a subclass
        // handles the list)
        const eris::eris_id_t author_ = 0;

        // The various comparison functions; one of these gets passed to std::stable_sort.
#define LESS_GREATER_METHODS(col) \
        static bool less_##col(const state::BookState &a, const state::BookState &b); \
        static bool greater_##col(const state::BookState &a, const state::BookState &b);
        LESS_GREATER_METHODS(id)
        LESS_GREATER_METHODS(author)
        LESS_GREATER_METHODS(pos_x)
        LESS_GREATER_METHODS(pos_y)
        LESS_GREATER_METHODS(quality)
        LESS_GREATER_METHODS(price)
        LESS_GREATER_METHODS(revenue)
        LESS_GREATER_METHODS(revenue_lifetime)
        LESS_GREATER_METHODS(pos_str)
        LESS_GREATER_METHODS(market_str)
        LESS_GREATER_METHODS(market_private)
        LESS_GREATER_METHODS(market_public)
        LESS_GREATER_METHODS(market_any)
        LESS_GREATER_METHODS(age)
        LESS_GREATER_METHODS(created)
        LESS_GREATER_METHODS(sales)
        LESS_GREATER_METHODS(sales_lifetime_private)
        LESS_GREATER_METHODS(sales_lifetime_public)
        LESS_GREATER_METHODS(sales_lifetime)
        LESS_GREATER_METHODS(pirated)
        LESS_GREATER_METHODS(pirated_lifetime)
        LESS_GREATER_METHODS(copies)
        LESS_GREATER_METHODS(copies_lifetime)
        LESS_GREATER_METHODS(lifetime_private)
#undef LESS_GREATER_METHODS
};

}}
