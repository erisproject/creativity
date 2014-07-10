#pragma once
#include <memory>
#include <gtkmm/treemodel.h>
#include <gtkmm/treesortable.h>
#include <gtkmm/treeview.h>
#include <gtkmm/treepath.h>
#include <eris/Simulation.hpp>

namespace creativity { namespace gui {

/** Base class for a flat list of simulation members.
 *
 * Note: the implementing class *must* inherit from Glib::ObjectBase and call the
 * Glib::ObjectBase(typeid(CLASS)) constructor to register the type with Glib.
 * */
template <class M>
class MemberStore : public Gtk::TreeModel, public Gtk::TreeSortable {
    public:
        MemberStore() = delete;

        /** Takes a Gtk::TreeView and adds this object's columns to it. */
        virtual void appendColumnsTo(Gtk::TreeView &v) const = 0;

        /** Gets the Member from a Path.  Returns an empty SharedMember if the Path is invalid. */
        eris::SharedMember<M> member(const Path &path) const;

        /** Resynchronizes the member list from the simulation.  This calls resync_remove() and
         * resync_changed, and resync_add() to get members to remove/add and notify of changes, then
         * takes the appropriate actions to remove/add them.  This also updates max_id_ if any of
         * the resync_add() members have an id larger than the current max_id_ value. */
        void resync();

    protected:
        /** Protected constructor; this object should be constructed from a subclass, typically via
         * the subclass's static create() method.
         *
         * \param sim a shared pointer value to the simulation object.
         *
         * \sa MemberStore
         * \sa BookStore
         */
        MemberStore(std::shared_ptr<eris::Simulation> &&sim);

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

        /// Returns a vector of new members to add to the model.  Called by resync().
        virtual std::vector<eris::SharedMember<M>> resync_add() = 0;

        /** Returns a set of members_ indices that should be removed from the model.  Called by
         * resync().  The default implementation returns any members that have an id() of 0 (which
         * indicates that a members doesn't belong to the simulation), but subclasses may override
         * this.
         */
        virtual std::unordered_set<size_t> resync_remove();

        /** Calls row_changed() on every row in the model.  Can be overridden to provide alternate
         * behaviour.
         */
        virtual void resync_changed();

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
                std::function<bool(const eris::SharedMember<M> &a, const eris::SharedMember<M> &b)> &compare,
                int sort_column_id,
                Gtk::SortType order);

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

        /** Can be used by subclasses to track the maximum member eris_id_t currently in members_,
         * so as to identify new members during resync().
         */
        eris::eris_id_t max_id_ = 0;
        /// The simulation.
        std::shared_ptr<eris::Simulation> sim_;
        /// The vector of members.
        std::vector<eris::SharedMember<M>> members_;
        /// Tracks model changes by being incremented whenever such a change occurs
        int stamp_ = 1;

    private:
        /// The sort column; the default is unsorted
        int sort_by_ = Gtk::TreeSortable::DEFAULT_UNSORTED_COLUMN_ID;
        /// The sort order; the default is ascending
        Gtk::SortType sort_order_ = Gtk::SORT_ASCENDING;
};

template <class M> MemberStore<M>::MemberStore(std::shared_ptr<eris::Simulation> &&sim)
    : //Glib::ObjectBase(typeid(MemberStore<M>)),
    //Gtk::TreeModel(),
    //Gtk::TreeSortable(),
    //Glib::Object(),
    sim_{std::move(sim)}
{}

template <class M> Gtk::TreeModelFlags MemberStore<M>::get_flags_vfunc() const {
    return Gtk::TREE_MODEL_LIST_ONLY;
}

template <class M> eris::SharedMember<M> MemberStore<M>::member(const Path &path) const {
    size_t i = path.back();
    if (i >= members_.size()) return eris::SharedMember<M>();
    return members_[i];
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
    size_t i = (size_t) iter.gobj()->user_data;
    if (i >= members_.size()) return false;

    iter_next.set_stamp(stamp_);
    iter_next.gobj()->user_data = (void*) ++i;
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
                std::function<bool(const eris::SharedMember<M> &a, const eris::SharedMember<M> &b)> &compare,
                int sort_column_id,
                Gtk::SortType order) {
    if (not compare) return;

    // The rows_reordered signal needs a vector where v[new_pos] == old_pos,
    // so store all the old positions before we sort.
    std::unordered_map<eris::eris_id_t, int> old_pos;
    int i = 0;
    for (auto &r : members_)
        old_pos[r] = i++;

    if (sort_by_ != sort_column_id or sort_order_ != order) {
        sort_by_ = sort_column_id;
        sort_order_ = order;
        sort_column_changed();
    }

    std::stable_sort(members_.begin(), members_.end(), compare);

    std::vector<int> new_order;
    new_order.reserve(members_.size());
    bool actual_reordering = false;
    for (auto &r : members_) {
        if (not actual_reordering and old_pos[r] != (int) new_order.size()) actual_reordering = true;
        new_order.push_back(old_pos[r]);
    }

    // If sorting actually changed anything, fire the rows-reordered signal
    if (actual_reordering) rows_reordered(Path(), new_order);
}

template <class M> void MemberStore<M>::resync() {
    auto remove = resync_remove();

    // If we're deleting, we have to copy out the whole vector, and will have to send out
    // row_deleted signals for the removed rows.
    if (not remove.empty()) {
        std::vector<eris::SharedMember<M>> pruned_readers;
        for (size_t i = 0; i < members_.size(); i++) {
            if (not remove.count(i)) {
                pruned_readers.push_back(std::move(members_[i]));
            }
        }
        // Swap the pruned list into the actual list and increment the validity stamp
        members_ = std::move(pruned_readers);
        stamp_++;

        // Send the row deleted signals
        for (auto &i : remove) {
            Path p;
            p.push_back(i);
            row_deleted(p);
        }
    }

    resync_changed();

    auto add = resync_add();

    if (not add.empty()) {
        size_t start_at = members_.size();
        members_.reserve(start_at + add.size());
        // Append the new readers
        for (auto &m : add) {
            if (m->id() > max_id_) max_id_ = m->id();
            members_.push_back(std::move(m));
        }

        if (remove.empty()) stamp_++; // If we removed some elements then we already incremented the stamp_

        // Send row_inserted signals for the new rows/readers
        for (size_t i = start_at; i < members_.size(); i++) {
            Path p;
            p.push_back(i);
            iterator it;
            it.gobj()->user_data = (void*) i;
            row_inserted(p, it);
        }

        // Finally resort the list (because the new elements won't be in the right locations)
        set_sort_column(sort_by_, sort_order_);
    }
}

template <class M> std::unordered_set<size_t> MemberStore<M>::resync_remove() {
    std::unordered_set<size_t> remove;
    size_t i = 0;
    for (auto &m : members_) {
        if (m == 0) remove.insert(i);
        i++;
    }
    return remove;
}

template <class M> void MemberStore<M>::resync_changed() {
    for (size_t i = 0; i < members_.size(); i++) {
        Path path;
        path.push_back(i);
        iterator iter;
        iter.set_stamp(stamp_);
        iter.gobj()->user_data = (void*) i;
        row_changed(path, iter);
    }
}

}}
