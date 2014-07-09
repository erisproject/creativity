#include "creativity/gui/ReaderStore.hpp"
#include "creativity/gui/GUI.hpp"
#include <gtkmm/treepath.h>
#include <algorithm>

using namespace eris;

namespace creativity { namespace gui {

    
ReaderStore::ReaderStore(std::shared_ptr<Simulation> &&sim)
    : Glib::ObjectBase(typeid(ReaderStore)),
    Gtk::TreeModel(),
    Gtk::TreeSortable(),
    Glib::Object(),
    sim_{std::move(sim)},
    readers_{sim_->agents<Reader>()}
{
    eris_id_t max_id = 0;
    for (auto &r : readers_)
        if (r->id() > max_id) max_id = r->id();

    max_id_ = max_id;
}

Glib::RefPtr<ReaderStore> ReaderStore::create(std::shared_ptr<Simulation> sim) {
    return Glib::RefPtr<ReaderStore>(new ReaderStore(std::move(sim)));
}

void ReaderStore::resync() {
    ERIS_DBG("resyncing");
    bool deleted = false;
    for (auto &r : readers_) {
        if (r == 0) {
            deleted = true;
            break;
        }
    }
    std::list<size_t> removed;
    // If we found any deleted, we have to copy out the whole vector, and will have to send out
    // row_deleted signals for the removed rows.
    if (deleted) {
        std::vector<SharedMember<Reader>> pruned_readers;
        for (size_t i = 0; i < readers_.size(); i++) {
            if (readers_[i] == 0) {
                removed.push_back(i);
            }
            else {
                pruned_readers.push_back(std::move(readers_[i]));
            }
        }
        // Swap the pruned list into the actual list and increment the validity stamp
        readers_.swap(pruned_readers);
        stamp_++;

        // Send the row deleted signals
        for (auto &i : removed) {
            Path p;
            p.push_back(i);
            row_deleted(p);
        }
    }

    // Look for new readers
    eris_id_t new_max = 0;
    auto new_readers = sim_->agents<Reader>([this,&new_max](const Reader &r) -> bool {
            if (r.id() > max_id_) {
                if (r.id() > new_max) new_max = r.id();
                return true;
            }
            return false;
    });
    if (new_max > 0) {
        ERIS_DBG("readers_ increasing from " << readers_.size() << " by " << new_readers.size() << " to " << readers_.size() + new_readers.size());
        size_t start_at = readers_.size();
        // Append the new readers
        readers_.insert(readers_.end(), new_readers.begin(), new_readers.end());
        max_id_ = new_max;

        if (not deleted) stamp_++; // If we also deleted, no need to increment again

        // Send row_inserted signals for the new rows/readers
        for (size_t i = start_at; i < readers_.size(); i++) {
            Path p;
            p.push_back(i);
            iterator it;
            it.gobj()->user_data = (void*) i;
            row_inserted(p, it);
        }
    }
}

Gtk::TreeModelFlags ReaderStore::get_flags_vfunc() const {
    return Gtk::TREE_MODEL_LIST_ONLY;
}

int ReaderStore::get_n_columns_vfunc() const {
    return columns.size();
}

GType ReaderStore::get_column_type_vfunc(int index) const {
    if ((unsigned int) index >= columns.size()) throw std::out_of_range("Invalid column index accessed");
    return columns.types()[index];
}

SharedMember<Reader> ReaderStore::reader(const Path &path) const {
    size_t i = path.back();
    if (i >= readers_.size()) return SharedMember<Reader>();
    return readers_[i];
}

bool ReaderStore::get_iter_vfunc(const Path &path, iterator& iter) const {
    size_t i = path.back();
    if (i >= readers_.size()) return false;

    iter.set_stamp(stamp_);
    iter.gobj()->user_data = (void*) i;
    return true;
}

bool ReaderStore::iter_next_vfunc(const iterator &iter, iterator &iter_next) const {
    if (iter.get_stamp() != stamp_) return false;
    size_t i = (size_t) iter.gobj()->user_data;
    if (i >= readers_.size()) return false;

    iter_next.set_stamp(stamp_);
    iter_next.gobj()->user_data = (void*) ++i;
    return true;
}

bool ReaderStore::iter_children_vfunc(const iterator&, iterator&) const {
    return false; // No nesting
}

bool ReaderStore::iter_parent_vfunc(const iterator&, iterator&) const {
    return false; // No nesting
}

bool ReaderStore::iter_nth_child_vfunc(const iterator&, int, iterator&) const {
    return false; // No nesting
}

bool ReaderStore::iter_nth_root_child_vfunc(int n, iterator &iter) const {
    size_t i = n;
    if (i >= readers_.size()) return false;

    iter.set_stamp(stamp_);
    iter.gobj()->user_data = (void*) i;
    return true;
}

bool ReaderStore::iter_has_child_vfunc(const iterator &) const {
    return false; // No nesting
}

int ReaderStore::iter_n_children_vfunc(const iterator &) const {
    return 0; // No nesting
}

int ReaderStore::iter_n_root_children_vfunc() const {
    return (int) readers_.size();
}

Gtk::TreeModel::Path ReaderStore::get_path_vfunc(const iterator &iter) const {
    size_t i = (size_t) iter.gobj()->user_data;
    Gtk::TreeModel::Path ret;
    ret.push_back(i);
    return ret;
}

void ReaderStore::get_value_vfunc(const iterator &iter, int column, Glib::ValueBase &value) const {
    if (iter.get_stamp() != stamp_ or column > get_n_columns_vfunc()) return;
    ERIS_DBG("(id=" << (size_t) iter.gobj()->user_data << ", col=" << column << ")");

    auto &r = readers_.at((size_t) iter.gobj()->user_data);
    if (column == columns.id.index()) {
        Glib::Value<eris_id_t> v;
        v.init(v.value_type());
        v.set(r->id());
        value.init(v.gobj());
    }
    else if (column == columns.posX.index() or column == columns.posY.index() or column == columns.u.index() or column == columns.uLifetime.index()) {
        Glib::Value<double> v;
        v.init(v.value_type());
        v.set(  column == columns.posX.index() ? r->position()[0] :
                column == columns.posY.index() ? r->position()[1] :
                column == columns.u.index() ? r->u() :
                r->uLifetime()
             );
        value.init(v.gobj());
    }
    else if (column == columns.posstr.index()) {
        Glib::Value<std::string> v;
        v.init(v.value_type());
        v.set(GUI::pos_to_string(r->position()));
        value.init(v.gobj());
    }
    else if (column == columns.booksOwned.index() or column == columns.booksNew.index() or column == columns.booksWritten.index()
            or column == columns.lastBookAge.index()) {
        Glib::Value<size_t> v;
        v.init(v.value_type());
        v.set(  column == columns.booksOwned.index() ? r->library().size() :
                column == columns.booksNew.index() ? r->newBooks().size() :
                column == columns.booksWritten.index() ? r->wrote().size() :
                r->wrote().empty() ? sim_->t() : r->wrote().back()->age()
             );
        value.init(v.gobj());
    }
    else {
        throw std::out_of_range("Invalid column index accessed");
    }
}

bool ReaderStore::get_sort_column_id_vfunc(int *sort_column_id, Gtk::SortType *order) const {
    if (sort_column_id) *sort_column_id = sort_by_;
    if (order) *order = sort_order_;
    return not(sort_by_ == DEFAULT_UNSORTED_COLUMN_ID or sort_by_ == DEFAULT_SORT_COLUMN_ID);
}

void ReaderStore::set_sort_column_id_vfunc(int sort_column_id, Gtk::SortType order) {
    ERIS_DBGVAR(sort_column_id);
    bool ascending = (order == Gtk::SORT_ASCENDING);
    ERIS_DBGVAR(ascending);
    std::function<bool(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b)> compare;
    if (sort_column_id == columns.id.index() || sort_column_id == DEFAULT_SORT_COLUMN_ID)
        compare = ascending ? less_id : greater_id;
    else if (sort_column_id == columns.posX.index())
        compare = ascending ? less_posX : greater_posX;
    else if (sort_column_id == columns.posY.index())
        compare = ascending ? less_posY : greater_posY;
    else if (sort_column_id == columns.posstr.index())
        compare = ascending ? less_posstr : greater_posstr;
    else if (sort_column_id == columns.u.index())
        compare = ascending ? less_uCurr : greater_uCurr;
    else if (sort_column_id == columns.uLifetime.index())
        compare = ascending ? less_uLife : greater_uLife;
    else if (sort_column_id == columns.booksOwned.index())
        compare = ascending ? less_booksOwned : greater_booksOwned;
    else if (sort_column_id == columns.booksNew.index())
        compare = ascending ? less_booksNew : greater_booksNew;
    else if (sort_column_id == columns.booksWritten.index())
        compare = ascending ? less_booksWritten : greater_booksWritten;
    else if (sort_column_id == columns.lastBookAge.index())
        compare = ascending ? less_lastBookAge : greater_lastBookAge;

    if (compare) {
        // The rows_reordered signal needs a vector where v[new_pos] == old_pos,
        // so store all the old positions before we sort.
        std::unordered_map<eris_id_t, int> old_pos;
        int i = 0;
        for (auto &r : readers_)
            old_pos[r] = i++;

        sort_by_ = sort_column_id;
        sort_order_ = order;
        sort_column_changed();

        std::stable_sort(readers_.begin(), readers_.end(), compare);

        std::vector<int> new_order;
        new_order.reserve(readers_.size());
        for (auto &r : readers_)
            new_order.push_back(old_pos[r]);

        rows_reordered(Path(), new_order);
    }
}

bool ReaderStore::less_id(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b) {
    return a->id() < b->id();
}
bool ReaderStore::greater_id(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b) {
    return a->id() > b->id();
}
bool ReaderStore::less_posX(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b) {
    return a->position()[0] < b->position()[0];
}
bool ReaderStore::greater_posX(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b) {
    return a->position()[0] > b->position()[0];
}
bool ReaderStore::less_posY(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b) {
    return a->position()[1] < b->position()[1];
}
bool ReaderStore::greater_posY(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b) {
    return a->position()[1] > b->position()[1];
}
// First x, then y for ties
bool ReaderStore::less_posstr(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b) {
    auto ax = a->position()[0], bx = b->position()[0];
    return ax == bx ? a->position()[1] < b->position()[1] : ax < bx;
}
bool ReaderStore::greater_posstr(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b) {
    auto ax = a->position()[0], bx = b->position()[0];
    return ax == bx ? a->position()[1] > b->position()[1] : ax > bx;
}
bool ReaderStore::less_uCurr(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b) {
    return a->u() < b->u();
}
bool ReaderStore::greater_uCurr(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b) {
    return a->u() > b->u();
}
bool ReaderStore::less_uLife(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b) {
    return a->uLifetime() < b->uLifetime();
}
bool ReaderStore::greater_uLife(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b) {
    return a->uLifetime() > b->uLifetime();
}
bool ReaderStore::less_booksOwned(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b) {
    return a->library().size() < b->library().size();
}
bool ReaderStore::greater_booksOwned(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b) {
    return a->library().size() > b->library().size();
}
bool ReaderStore::less_booksNew(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b) {
    return a->newBooks().size() < b->newBooks().size();
}
bool ReaderStore::greater_booksNew(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b) {
    return a->newBooks().size() > b->newBooks().size();
}
bool ReaderStore::less_booksWritten(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b) {
    return a->wrote().size() < b->wrote().size();
}
bool ReaderStore::greater_booksWritten(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b) {
    return a->wrote().size() > b->wrote().size();
}
bool ReaderStore::less_lastBookAge(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b) {
    return (a->wrote().empty() ? std::numeric_limits<unsigned long>::max() : a->wrote().back()->age())
         < (b->wrote().empty() ? std::numeric_limits<unsigned long>::max() : b->wrote().back()->age());
}
bool ReaderStore::greater_lastBookAge(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b) {
    return (a->wrote().empty() ? std::numeric_limits<unsigned long>::max() : a->wrote().back()->age())
         > (b->wrote().empty() ? std::numeric_limits<unsigned long>::max() : b->wrote().back()->age());
}

void ReaderStore::appendColumnsTo(Gtk::TreeView &v) const {
    appendCol(v, "ID", columns.id, 100);
    appendCol(v, "Position", columns.posstr, 150);
    appendCol(v, "Utility", columns.u, 100);
    appendCol(v, "Life Util.", columns.uLifetime, 100);
    appendCol(v, "Books", columns.booksOwned, 100);
    appendCol(v, "# New", columns.booksNew, 100);
    appendCol(v, "# Written", columns.booksWritten, 100);
    appendCol(v, "Last wrote", columns.lastBookAge, 100);
}

}}
