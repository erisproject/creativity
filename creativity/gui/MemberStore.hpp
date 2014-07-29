#pragma once
#include <memory>
#include <gtkmm/treemodel.h>
#include <gtkmm/treesortable.h>
#include <gtkmm/treeview.h>
#include <gtkmm/treepath.h>
#include "creativity/state/State.hpp"

namespace creativity { namespace gui {

/** Base class for a flat list of simulation members.
 *
 * Note: the implementing class *must* inherit from Glib::ObjectBase and call the
 * Glib::ObjectBase(typeid(CLASS)) constructor to register the type with Glib.
 *
 * \param M the state type being represented by this MemberStore such as `BookState` or
 * `ReaderState`.
 */
template <class M>
class MemberStore : public Gtk::TreeModel, public Gtk::TreeSortable {
    public:
        MemberStore() = delete;

        /** Takes a Gtk::TreeView and adds this object's columns to it. */
        virtual void appendColumnsTo(Gtk::TreeView &v) const = 0;

        /** Gets the M Member from a Path.  Throws an exception if the Path is invalid. */
        const M& member(const Path &path) const;

        /** Gets the M Member from an iterator.  Throws an exception if the iterator is invalid. */
        const M& member(const iterator &iter) const;

        /** Returns the Path to the member with the given id.  Returns an empty path if the id was
         * not found.
         *
         * This does a linear search of the current member list.  If the optional `hint` is given,
         * locations near `hint` are checked first; otherwise this does a linear search beginning
         * from 0.
         *
         * \param id the id to compare against members_[i].id to find the element's position.
         * \param hint if provided as a non-zero, this starts looking at and near the given
         * position, iterating outward in both directions.  This is a useful optimization when
         * calling find() on a model that has a reasonable chance of having the elements in the same
         * (or nearby) positions.
         */
        Path find(eris::eris_id_t id, size_t hint = 0) const;

        /** Returns the Path to the member with the given id.  The `hint` iterator is used to
         * determine the position to check first; if it doesn't match the requested ID, a linear
         * search of all members is performed.
         *
         * \sa find(eris::eris_id_t, size_t)
         */
        Path find(eris::eris_id_t id, const iterator &iter) const;

        /** Returns the Path to the member with the given id.  The `hint` path is used to determine
         * the position to check first; if it doesn't match the requested ID, a linear search of all
         * members is performed.
         *
         * \sa find(eris::eris_id_t, size_t)
         */
        Path find(eris::eris_id_t id, const Path &hint) const;

    protected:
        /** Protected constructor; this object should be constructed from a subclass, typically via
         * the subclass's static create() method.
         *
         * Subclasses must inherit from Glib::ObjectBase and call the
         * Glib::ObjectBase(typeid(CLASS)) constructor to register the type with Glib.  See
         * gui/ReaderStore.hpp for an example.
         *
         * \param state a reference to a simulation state object.  The object's lifespan should
         * exceed the lifespan of the MemberStore-derived object.
         *
         * \sa ReaderStore
         * \sa BookStore
         */
        MemberStore(std::shared_ptr<const state::State> &&state);

        /** Returns Gtk::TreeModel flags (specifically, the LIST_ONLY flag). */
        virtual Gtk::TreeModelFlags get_flags_vfunc() const override;

        /** Returns `obj.columns.size()`, the number of model columns. */
        virtual int get_n_columns_vfunc() const override = 0;

        /** Returns the column type of the given position.  This is typically invoked via
         * get_column_type, itself given a column member of the `.columns` ColRec object.
         *
         * \sa MemberStore::ColRec
         * \sa BookStore::ColRec
         */
        virtual GType get_column_type_vfunc(int index) const override = 0;

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
        /// Returns false always: MemberStore elements cannot have children.
        virtual bool iter_children_vfunc(const iterator&, iterator&) const override;
        /// Returns false always: MemberStore elemenets cannot have children/parents.
        virtual bool iter_parent_vfunc(const iterator&, iterator&) const override;
        /// Returns false always: MemberStore elements cannot have children.
        virtual bool iter_nth_child_vfunc(const iterator&, int, iterator&) const override;
        /// Returns false always: MemberStore elements cannot have children.
        virtual bool iter_has_child_vfunc(const iterator &) const override;
        /// Returns 0 always: MemberStore elements have no children.
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
        virtual void get_value_vfunc(const iterator &iter, int column, Glib::ValueBase &value) const override = 0;

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
         */
        virtual void set_sort_column_id_vfunc(int sort_column_id, Gtk::SortType order) override = 0;

        /** Called by subclasses, typically in set_sort_column_id_vfunc, to resort the current list
         * of members using the given function.  This uses std::stable_sort, which means the current
         * order is maintained for any members that are equal according to the given comparison
         * function.
         *
         * This fires sort_column_changed if the sort column and/or sort order are changing from
         * their current values, and rows_reordered if the sorting actually changes any row orders.
         *
         * If compare is not callable, this method does nothing.
         *
         * \param compare a function (or callable object) that returns true if its first arguments
         * should (strictly) come before its second argument, and false otherwise.
         * \param sort_column_id the new sort column
         * \param order the new sort order
         */
        virtual void sort_members(
                std::function<bool(const M &a, const M &b)> &compare,
                int sort_column_id,
                Gtk::SortType order);

        /** Appends a single column to the given view using the given label, width, and sortability.
         *
         * \param v the TreeView to add the column to
         * \param label the string for the column header
         * \param col the model column object, typically a member of `store.columns` such as
         * `readerstore.columns.id`.
         * \param width the fixed width of the column
         * \param sortable if omitted (or true) the column will be user-sortable by clicking on the
         * column header; if specified explicitly as false, user sorting will be disabled.
         */
        template <typename T, typename = typename std::enable_if<std::is_base_of<Gtk::TreeModelColumnBase, T>::value>>
        void appendCol(Gtk::TreeView &v, const std::string &label, T &col, int width, bool sortable = true) const {
            v.append_column(label, col);
            auto *c = v.get_column(v.get_n_columns()-1);
            c->set_sizing(Gtk::TREE_VIEW_COLUMN_FIXED);
            c->set_fixed_width(width);
            if (sortable)
                c->set_sort_column(col);
        }

        /// The state this MemberStore represents.
        const std::shared_ptr<const state::State> state_;
        /// The vector of members.  Subclasses need to add all member references to this vector.
        std::vector<std::reference_wrapper<const M>> members_;
        /// Tracks model changes by being incremented whenever such a change occurs
        int stamp_ = 1;

    private:
        /// The sort column; the default is unsorted
        int sort_by_ = Gtk::TreeSortable::DEFAULT_UNSORTED_COLUMN_ID;
        /// The sort order; the default is ascending
        Gtk::SortType sort_order_ = Gtk::SORT_ASCENDING;
};

template <class M> MemberStore<M>::MemberStore(std::shared_ptr<const state::State> &&state) : state_{std::move(state)} {}

template <class M> Gtk::TreeModelFlags MemberStore<M>::get_flags_vfunc() const {
    return Gtk::TREE_MODEL_LIST_ONLY;
}

template <class M> const M& MemberStore<M>::member(const Path &path) const {
    return members_.at(path.back());
}

template <class M> const M& MemberStore<M>::member(const iterator &iter) const {
    return members_.at((size_t) iter.gobj()->user_data);
}

template <class M> Gtk::TreeModel::Path MemberStore<M>::find(eris::eris_id_t id, size_t hint) const {
    Gtk::TreeModel::Path p;
    const long max = members_.size();
    // Loop through members_ starting at hint with two counters: one starting at hint and going down
    // to 0, the other starting at hint+1 and going up to max.  This makes us search positions
    // closer to hint first.
    long a = hint, b = hint + 1;
    while (b < max or a >= 0) {
        if (a >= 0) {
            if (members_[a].get().id == id) {
                p.push_back(a);
                break;
            }
            a--;
        }
        if (b < max) {
            if (members_[b].get().id == id) {
                p.push_back(b);
                break;
            }
            b++;
        }
    }
    return p;
}
template <class M> Gtk::TreeModel::Path MemberStore<M>::find(eris::eris_id_t id, const iterator &hint) const {
    return find(id, (size_t) hint.gobj()->user_data);
}
template <class M> Gtk::TreeModel::Path MemberStore<M>::find(eris::eris_id_t id, const Path &hint) const {
    return find(id, hint.empty() ? 0 : hint.back());
}

template <class M> bool MemberStore<M>::get_iter_vfunc(const Path &path, iterator& iter) const {
    size_t i = path.back();
    if (i >= members_.size()) return false;

    iter.set_stamp(stamp_);
    iter.gobj()->user_data = (void*) i;
    return true;
}

template <class M> bool MemberStore<M>::iter_next_vfunc(const iterator &iter, iterator &iter_next) const {
    if (iter.get_stamp() != stamp_) return false;
    size_t i = 1 + (size_t) iter.gobj()->user_data;
    if (i >= members_.size()) return false;

    iter_next.set_stamp(stamp_);
    iter_next.gobj()->user_data = (void*) i;
    return true;
}

template <class M> bool MemberStore<M>::iter_children_vfunc(const iterator&, iterator&) const {
    return false; // No nesting
}

template <class M> bool MemberStore<M>::iter_parent_vfunc(const iterator&, iterator&) const {
    return false; // No nesting
}

template <class M> bool MemberStore<M>::iter_nth_child_vfunc(const iterator&, int, iterator&) const {
    return false; // No nesting
}

template <class M> bool MemberStore<M>::iter_nth_root_child_vfunc(int n, iterator &iter) const {
    size_t i = n;
    if (i >= members_.size()) return false;

    iter.set_stamp(stamp_);
    iter.gobj()->user_data = (void*) i;
    return true;
}

template <class M> bool MemberStore<M>::iter_has_child_vfunc(const iterator &) const {
    return false; // No nesting
}

template <class M> int MemberStore<M>::iter_n_children_vfunc(const iterator &) const {
    return 0; // No nesting
}

template <class M> int MemberStore<M>::iter_n_root_children_vfunc() const {
    return (int) members_.size();
}

template <class M> Gtk::TreeModel::Path MemberStore<M>::get_path_vfunc(const iterator &iter) const {
    size_t i = (size_t) iter.gobj()->user_data;
    Gtk::TreeModel::Path ret;
    ret.push_back(i);
    return ret;
}

template <class M> bool MemberStore<M>::get_sort_column_id_vfunc(int *sort_column_id, Gtk::SortType *order) const {
    if (sort_column_id) *sort_column_id = sort_by_;
    if (order) *order = sort_order_;
    return not(sort_by_ == DEFAULT_UNSORTED_COLUMN_ID or sort_by_ == DEFAULT_SORT_COLUMN_ID);
}

template <class M> void MemberStore<M>::sort_members(
                std::function<bool(const M &a, const M &b)> &compare,
                int sort_column_id,
                Gtk::SortType order) {
    if (not compare) return;

    // The rows_reordered signal needs a vector where v[new_pos] == old_pos,
    // so store all the old positions before we sort.
    std::unordered_map<eris::eris_id_t, int> old_pos;
    int i = 0;
    for (const M &m : members_)
        old_pos[m.id] = i++;

    if (sort_by_ != sort_column_id or sort_order_ != order) {
        sort_by_ = sort_column_id;
        sort_order_ = order;
        sort_column_changed();
    }

    std::stable_sort(members_.begin(), members_.end(), compare);

    // The sort invalidates our tree iterators, so increment the stamp
    stamp_++;

    std::vector<int> new_order;
    new_order.reserve(members_.size());
    bool actual_reordering = false;
    for (const M &m : members_) {
        if (not actual_reordering and old_pos[m.id] != (int) new_order.size()) actual_reordering = true;
        new_order.push_back(old_pos[m.id]);
    }

    // If sorting actually changed anything, fire the rows-reordered signal
    if (actual_reordering) rows_reordered(Path(), new_order);
}

}}
