#include "creativity/gui/ReaderStore.hpp"
#include "creativity/gui/GUI.hpp"
#include <gtkmm/treepath.h>

using namespace eris;

namespace creativity { namespace gui {

    
ReaderStore::ReaderStore(std::shared_ptr<Simulation> &&sim)
    : Glib::ObjectBase(typeid(ReaderStore)),
    Gtk::TreeModel(),
    Glib::Object(),
    sim_{sim},
    readers_{sim->agents<Reader>()}
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

// void ReaderStore::ref_node_vfunc(const iterator &iter) const {}
//void ReaderStore::unref_node_vfunc(const iterator &iter) const {}

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
            or column == columns.bookLastAge.index()) {
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

void ReaderStore::appendColumnsTo(Gtk::TreeView &v) const {
    appendCol(v, "ID", columns.id, 100);
    appendCol(v, "Position", columns.posstr, 150, false);
    appendCol(v, "Utility", columns.u, 100);
    appendCol(v, "Life Util.", columns.uLifetime, 100);
    appendCol(v, "Books", columns.booksOwned, 100);
    appendCol(v, "# New", columns.booksNew, 100);
    appendCol(v, "# Written", columns.booksWritten, 100);
    appendCol(v, "Last wrote", columns.bookLastAge, 100);
}

}}
