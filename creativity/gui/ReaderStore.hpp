#pragma once
#include <memory>
#include <gtkmm/treemodel.h>
#include <gtkmm/treesortable.h>
#include <gtkmm/treeview.h>
#include <eris/Simulation.hpp>
#include "creativity/Reader.hpp"

namespace creativity { namespace gui {

/** Gtk::TreeModel::ColumnRecord subclass for handling Reader information in the list of readers in
 * the main GUI window.
 */
class ReaderStore : public Gtk::TreeModel, public Gtk::TreeSortable, Glib::Object {
    public:
        ReaderStore() = delete;

        /** Interface class between a simulation's Readers and a Gtk::TreeView.  This internally
         * stores a vector of Readers which can be updated (if needed) by calling the update method.
         *
         * This exposes 9 columns:
         * - ID
         * - x position
         * - y position
         * - position string (made from x position and y position)
         * - current utility
         * - lifetime utility
         * - books owned
         * - new books
         * - books written
         * - age of most recently written book (simulation age if no books written)
         */
        static Glib::RefPtr<ReaderStore> create(std::shared_ptr<eris::Simulation> sim);

        /** Sychronizes the list of readers with the stored Simulation. Typically called after a
         * simulation period runs. */
        void resync();

        /** ColumnRecord object for a ReaderStore.  This object contains the columns for this Book
         * model.  This should not be used directly, but rather accessed via the public `columns`
         * member.
         */
        class ColRec : public Gtk::TreeModel::ColumnRecord {
            public:
                Gtk::TreeModelColumn<eris::eris_id_t> id;
                Gtk::TreeModelColumn<double> posX, posY, u, uLifetime;
                Gtk::TreeModelColumn<std::string> posstr;
                Gtk::TreeModelColumn<size_t> booksOwned, booksNew, booksWritten, lastBookAge;

            private:
                ColRec() {
                    add(id); add(posX); add(posY); add(posstr); add(u); add(uLifetime);
                    add(booksOwned); add(booksNew); add(booksWritten); add(lastBookAge);
                }
                friend class ReaderStore;
        };

        /** The columns of this ReaderStore.  For example, to access the utility column, use
         * `readerstore.columns.u`.
         */
        ColRec columns;

        /** Takes a Gtk::TreeView and adds this object's columns to it. */
        void appendColumnsTo(Gtk::TreeView &v) const;

        /** Gets a Reader from a Path.  Returns an empty SharedMember if the Path is invalid. */
        eris::SharedMember<Reader> reader(const Path &path) const;

    protected:
        /// Protected constructor; this object should be constructed using create().
        ReaderStore(std::shared_ptr<eris::Simulation> &&sim);

        /** Returns Gtk::TreeModel flags (specifically, the LIST_ONLY flag). */
        virtual Gtk::TreeModelFlags get_flags_vfunc() const override;
        /** Returns `obj.columns.size()`, the number of reader model columns. */
        virtual int get_n_columns_vfunc() const override;
        /** Returns the column type of the given position.  See the list of virtual columns in the
         * class description.
         *
         * \sa ReaderStore::ColRec
         */
        virtual GType get_column_type_vfunc(int index) const override;

        /** Converts a path to an iterator.
         *
         * \param path the path (in)
         * \param iter the iterator to set (out)
         * \returns true and sets `iter` if the path refers to a valid element, returns false otherwise.
         */
        virtual bool get_iter_vfunc(const Path &path, iterator& iter) const override;
        /** Takes an iterator, returns an iterator to the next item.
         *
         * \param iter the current element (in)
         * \param iter_next the next element (out)
         * \returns true and sets `iter_next` if `iter` is valid (i.e. the model hasn't changed
         * since the iterator was created) and there is a next element; returns false otherwise.
         */
        virtual bool iter_next_vfunc(const iterator &iter, iterator &iter_next) const override;
        /// Returns false always: ReaderStore elements cannot have children.
        virtual bool iter_children_vfunc(const iterator&, iterator&) const override;
        /// Returns false always: ReaderStore elemenets cannot have children/parents.
        virtual bool iter_parent_vfunc(const iterator&, iterator&) const override;
        /// Returns false always: ReaderStore elements cannot have children.
        virtual bool iter_nth_child_vfunc(const iterator&, int, iterator&) const override;
        /// Returns false always: ReaderStore elements cannot have children.
        virtual bool iter_has_child_vfunc(const iterator &) const override;
        /// Returns 0 always: ReaderStore elements have no children.
        virtual int iter_n_children_vfunc(const iterator &) const override;
        /** Obtains an iterator to the `n`th reader.
         *
         * \param n the index of the reader to access
         * \param iter an iterator to set to the requested reader
         * \returns true and sets iter if `n` is valid (i.e. there are at least `n+1` readers);
         * false otherwise.
         */
        virtual bool iter_nth_root_child_vfunc(int n, iterator &iter) const override;
        /// Returns the number of readers stored in this model.
        virtual int iter_n_root_children_vfunc() const override;
        /// Converts iterator `iter` into a Path.
        virtual Path get_path_vfunc(const iterator &iter) const override;
        /** Accesses a column value.
         *
         * \param iter a valid iterator referencing the row to access
         * \param column the index of the column to access
         * \param value a Glib::Value<TYPE> object (where TYPE is the appropriate type for the
         * requested `column`) in which to store the value.
         */
        virtual void get_value_vfunc(const iterator &iter, int column, Glib::ValueBase &value) const override;

        // TreeSortable overrides:

        /** Accesses the current sort column and order.  Returns true if accessed sort_column_id
         * refers to a specific column (instead of the Gtk magic unsorted and default order
         * constants).
         *
         * \param sort_column_id a pointer in which to store the sort column.  May be a nullptr to
         * skip accessing the sort column.
         * \param order a pointer in which to store the sort order.  May be a nullptr to skip
         * accessing the sort order.
         */
        virtual bool get_sort_column_id_vfunc(int *sort_column_id, Gtk::SortType *order) const override;

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
        // The maximum reader eris_id_t currently in readers_, used to identify new readers during resync():
        eris::eris_id_t max_id_;
        std::shared_ptr<eris::Simulation> sim_;
        std::vector<eris::SharedMember<Reader>> readers_;
        // Tracks model changes by being incremented whenever such a change occurs
        int stamp_ = 1;
        // Default is unsorted
        int sort_by_ = Gtk::TreeSortable::DEFAULT_UNSORTED_COLUMN_ID;
        Gtk::SortType sort_order_ = Gtk::SORT_ASCENDING;

        // The various comparison functions; one of these gets passed to std::stable_sort.
        static bool less_id(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool greater_id(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool less_posX(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool greater_posX(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool less_posY(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool greater_posY(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool less_posstr(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool greater_posstr(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool less_uCurr(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool greater_uCurr(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool less_uLife(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool greater_uLife(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool less_booksOwned(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool greater_booksOwned(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool less_booksNew(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool greater_booksNew(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool less_booksWritten(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool greater_booksWritten(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool less_lastBookAge(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool greater_lastBookAge(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);

        // Appends a single column to the given view using the given label, width, and sortability.
        template <typename T, typename = typename std::enable_if<std::is_base_of<Gtk::TreeModelColumnBase, T>::value>>
        void appendCol(Gtk::TreeView &v, const std::string &label, T &col, int width, bool sortable = true) const {
            v.append_column(label, col);
            auto *c = v.get_column(v.get_n_columns()-1);
            c->set_sizing(Gtk::TREE_VIEW_COLUMN_FIXED);
            c->set_fixed_width(width);
            if (sortable)
                c->set_sort_column(col);
        }

};

}}
