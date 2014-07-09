#include "creativity/gui/BookStore.hpp"
#include "creativity/gui/GUI.hpp"
#include <gtkmm/treepath.h>
#include <algorithm>

using namespace eris;
using namespace std::placeholders;

namespace creativity { namespace gui {
    
BookStore::BookStore(std::shared_ptr<Simulation> &&sim, SharedMember<Reader> &&author)
    : Glib::ObjectBase(typeid(BookStore)),
    Gtk::TreeModel(),
    Gtk::TreeSortable(),
    Glib::Object(),
    sim_{std::move(sim)},
    author_{std::move(author)}
{
    eris_id_t max_id = 0;
    if (author_.ptr()) {
        // Author list
        books_.reserve(author_->wrote().size());
        for (auto &b : author_->wrote()) {
            books_.push_back(b);
            if (b > max_id) max_id = b;
        }
    }
    else {
        // Simulation list
        books_ = sim_->goods<Book>();
        for (auto &b : books_)
            if (b > max_id) max_id = b;
    }

    max_id_ = max_id;
}

Glib::RefPtr<BookStore> BookStore::create(std::shared_ptr<Simulation> sim, SharedMember<Reader> author) {
    return Glib::RefPtr<BookStore>(new BookStore(std::move(sim), std::move(author)));
}

void BookStore::resync() {
    ERIS_DBG("resyncing");
    copies_cache_.clear();

    bool deleted = false;
    for (auto &b : books_) {
        if (b == 0) { // Removed from simulation
            deleted = true;
            break;
        }
    }
    std::list<size_t> removed;
    // If we found any deleted, we have to copy out the whole vector, and will have to send out
    // row_deleted signals for the removed rows.
    if (deleted) {
        std::vector<SharedMember<Book>> pruned_books;
        for (size_t i = 0; i < books_.size(); i++) {
            if (books_[i] == 0) {
                removed.push_back(i);
            }
            else {
                pruned_books.push_back(std::move(books_[i]));
            }
        }
        // Swap the pruned list into the actual list and increment the validity stamp
        books_.swap(pruned_books);
        stamp_++;

        // Send the row deleted signals
        for (auto &i : removed) {
            Path p;
            p.push_back(i);
            row_deleted(p);
        }
    }

    // Look for new books
    eris_id_t new_max = 0;
    std::vector<SharedMember<Book>> new_books;
    // NB: check .ptr because the author might have been removed from the simulation, but we are
    // still responsible for that author's books and don't want to revert to all simulation books.
    if (author_.ptr()) {
        for (auto &b : author_->wrote()) {
            if (b > max_id_) {
                if (b > new_max) new_max = b;
                new_books.push_back(b);
            }
        }
    }
    else {
        new_books = sim_->goods<Book>([this,&new_max](const Book &b) -> bool {
            if (b.id() > max_id_) {
                if (b.id() > new_max) new_max = b.id();
                return true;
            }
            return false;
        });
    }
    if (not new_books.empty()) {
        ERIS_DBG("books_ increasing from " << books_.size() << " by " << new_books.size() << " to " << books_.size() + new_books.size());
        size_t start_at = books_.size();
        // Append the new books
        books_.insert(books_.end(), new_books.begin(), new_books.end());
        max_id_ = new_max;

        if (not deleted) stamp_++; // If we also deleted, no need to increment again

        // Send row_inserted signals for the new rows/books
        for (size_t i = start_at; i < books_.size(); i++) {
            Path p;
            p.push_back(i);
            iterator it;
            it.gobj()->user_data = (void*) i;
            row_inserted(p, it);
        }
    }
}

Gtk::TreeModelFlags BookStore::get_flags_vfunc() const {
    return Gtk::TREE_MODEL_LIST_ONLY;
}

int BookStore::get_n_columns_vfunc() const {
    return columns.size();
}

GType BookStore::get_column_type_vfunc(int index) const {
    if ((unsigned int) index >= columns.size()) throw std::out_of_range("Invalid column index accessed");
    return columns.types()[index];
}

SharedMember<Book> BookStore::book(const Path &path) const {
    size_t i = path.back();
    if (i >= books_.size()) return SharedMember<Book>();
    return books_[i];
}

bool BookStore::get_iter_vfunc(const Path &path, iterator& iter) const {
    size_t i = path.back();
    if (i >= books_.size()) return false;

    iter.set_stamp(stamp_);
    iter.gobj()->user_data = (void*) i;
    return true;
}

bool BookStore::iter_next_vfunc(const iterator &iter, iterator &iter_next) const {
    if (iter.get_stamp() != stamp_) return false;
    size_t i = (size_t) iter.gobj()->user_data;
    if (i >= books_.size()) return false;

    iter_next.set_stamp(stamp_);
    iter_next.gobj()->user_data = (void*) ++i;
    return true;
}

bool BookStore::iter_children_vfunc(const iterator&, iterator&) const {
    return false; // No nesting
}

bool BookStore::iter_parent_vfunc(const iterator&, iterator&) const {
    return false; // No nesting
}

bool BookStore::iter_nth_child_vfunc(const iterator&, int, iterator&) const {
    return false; // No nesting
}

bool BookStore::iter_nth_root_child_vfunc(int n, iterator &iter) const {
    size_t i = n;
    if (i >= books_.size()) return false;

    iter.set_stamp(stamp_);
    iter.gobj()->user_data = (void*) i;
    return true;
}

bool BookStore::iter_has_child_vfunc(const iterator &) const {
    return false; // No nesting
}

int BookStore::iter_n_children_vfunc(const iterator &) const {
    return 0; // No nesting
}

int BookStore::iter_n_root_children_vfunc() const {
    return (int) books_.size();
}

Gtk::TreeModel::Path BookStore::get_path_vfunc(const iterator &iter) const {
    size_t i = (size_t) iter.gobj()->user_data;
    Gtk::TreeModel::Path ret;
    ret.push_back(i);
    return ret;
}

void BookStore::get_value_vfunc(const iterator &iter, int column, Glib::ValueBase &value) const {
    if (iter.get_stamp() != stamp_ or column > get_n_columns_vfunc()) return;
    ERIS_DBG("(id=" << (size_t) iter.gobj()->user_data << ", col=" << column << ")");

    auto &b = books_.at((size_t) iter.gobj()->user_data);
    if (column == columns.id.index() || column == columns.author.index() || column == columns.market.index()) {
        Glib::Value<eris_id_t> v;
        v.init(v.value_type());
        v.set(b->id());
        value.init(v.gobj());
    }
    else if (column == columns.posX.index() or column == columns.posY.index() or column == columns.quality.index() or column == columns.price.index()
            or column == columns.revenue.index() or column == columns.revenueLifetime.index()) {
        Glib::Value<double> v;
        v.init(v.value_type());
        v.set(  column == columns.posX.index() ? b->position()[0] :
                column == columns.posY.index() ? b->position()[1] :
                column == columns.quality.index() ? b->quality() :
                column == columns.price.index() ? (b->hasMarket() ? b->market()->price() : std::numeric_limits<double>::quiet_NaN()) :
                column == columns.revenue.index() ? b->currRevenue() :
                b->lifeRevenue()
             );
        value.init(v.gobj());
    }
    else if (column == columns.hasMarket.index()) {
        Glib::Value<gboolean> v;
        v.init(v.value_type());
        v.set(b->hasMarket());
        value.init(v.gobj());
    }
    else if (column == columns.posstr.index()) {
        Glib::Value<std::string> v;
        v.init(v.value_type());
        v.set(GUI::pos_to_string(b->position()));
        value.init(v.gobj());
    }
    else if (column == columns.copies.index()) {
        Glib::Value<size_t> v;
        v.init(v.value_type());
        if (copies_cache_.count(b) == 0) copies_cache_[b] = b->copies();
        v.set(copies_cache_[b]);
    }
    else if (column == columns.age.index() or column == columns.sales.index() or column == columns.salesLifetime.index()
            or column == columns.copies.index() or column == columns.lifetime.index()) {
        Glib::Value<size_t> v;
        v.init(v.value_type());
        v.set(  column == columns.age.index() ? b->age() :
                column == columns.sales.index() ? b->currSales() :
                column == columns.salesLifetime.index() ? b->lifeSales() :
                b->marketPeriods()
             );
        value.init(v.gobj());
    }
    else {
        throw std::out_of_range("Invalid column index accessed");
    }
}

bool BookStore::get_sort_column_id_vfunc(int *sort_column_id, Gtk::SortType *order) const {
    if (sort_column_id) *sort_column_id = sort_by_;
    if (order) *order = sort_order_;
    return not(sort_by_ == DEFAULT_UNSORTED_COLUMN_ID or sort_by_ == DEFAULT_SORT_COLUMN_ID);
}

void BookStore::set_sort_column_id_vfunc(int sort_column_id, Gtk::SortType order) {
    ERIS_DBGVAR(sort_column_id);
    bool ascending = (order == Gtk::SORT_ASCENDING);
    ERIS_DBGVAR(ascending);
    std::function<bool(const eris::SharedMember<Book> &a, const eris::SharedMember<Book> &b)> compare;
#define ELSE_IF_COL(COL) \
    else if (sort_column_id == columns.COL.index()) compare = ascending ? less_##COL : greater_##COL
    if (sort_column_id == columns.id.index() || sort_column_id == DEFAULT_SORT_COLUMN_ID)
        compare = ascending ? less_id : greater_id;
    ELSE_IF_COL(author);
    ELSE_IF_COL(market);
    ELSE_IF_COL(hasMarket);
    ELSE_IF_COL(posX);
    ELSE_IF_COL(posY);
    ELSE_IF_COL(posstr);
    ELSE_IF_COL(quality);
    ELSE_IF_COL(price);
    ELSE_IF_COL(revenue);
    ELSE_IF_COL(revenueLifetime);
    ELSE_IF_COL(age);
    ELSE_IF_COL(sales);
    ELSE_IF_COL(salesLifetime);
    ELSE_IF_COL(lifetime);
#undef ELSE_IF_COL
    else if (sort_column_id == columns.copies.index()) {
        // Can't handle like the above: less/greater_copies() is not static
        compare = std::bind(ascending ? &BookStore::less_copies : &BookStore::greater_copies, this, _1, _2);
    }

    if (compare) {
        // The rows_reordered signal needs a vector where v[new_pos] == old_pos,
        // so store all the old positions before we sort.
        std::unordered_map<eris_id_t, int> old_pos;
        int i = 0;
        for (auto &r : books_)
            old_pos[r] = i++;

        sort_by_ = sort_column_id;
        sort_order_ = order;
        sort_column_changed();

        std::stable_sort(books_.begin(), books_.end(), compare);

        std::vector<int> new_order;
        new_order.reserve(books_.size());
        for (auto &r : books_)
            new_order.push_back(old_pos[r]);

        rows_reordered(Path(), new_order);
    }
}

#define LESS_GREATER(COL, ACCESS) \
bool BookStore::less_##COL(const eris::SharedMember<Book> &a, const eris::SharedMember<Book> &b) { \
    return a->ACCESS < b->ACCESS; \
} \
bool BookStore::greater_##COL(const eris::SharedMember<Book> &a, const eris::SharedMember<Book> &b) { \
    return a->ACCESS > b->ACCESS; \
}
LESS_GREATER(id, id())
LESS_GREATER(posX, position()[0])
LESS_GREATER(posY, position()[1])
LESS_GREATER(author, author()->id())
LESS_GREATER(quality, quality())
LESS_GREATER(revenue, currRevenue())
LESS_GREATER(revenueLifetime, lifeRevenue())
LESS_GREATER(age, age())
LESS_GREATER(sales, currSales())
LESS_GREATER(salesLifetime, lifeSales())
LESS_GREATER(lifetime, marketPeriods())
LESS_GREATER(hasMarket, hasMarket())
#undef LESS_GREATER
// The rest are a bit more complicated:
bool BookStore::less_market(const eris::SharedMember<Book> &a, const eris::SharedMember<Book> &b) {
    return (a->hasMarket() ? a->market()->id() : 0) < (b->hasMarket() ? b->market()->id() : 0);
}
bool BookStore::greater_market(const eris::SharedMember<Book> &a, const eris::SharedMember<Book> &b) {
    return (a->hasMarket() ? a->market()->id() : 0) > (b->hasMarket() ? b->market()->id() : 0);
}
// First x, then y for ties
bool BookStore::less_posstr(const eris::SharedMember<Book> &a, const eris::SharedMember<Book> &b) {
    auto ax = a->position()[0], bx = b->position()[0];
    return ax == bx ? a->position()[1] < b->position()[1] : ax < bx;
}
bool BookStore::greater_posstr(const eris::SharedMember<Book> &a, const eris::SharedMember<Book> &b) {
    auto ax = a->position()[0], bx = b->position()[0];
    return ax == bx ? a->position()[1] > b->position()[1] : ax > bx;
}
// For sorting purposes, marketless books are considered to have a price of negative infinity
// (rather than the model value of quiet_NaN)
bool BookStore::less_price(const eris::SharedMember<Book> &a, const eris::SharedMember<Book> &b) {
    return (a->hasMarket() ? a->market()->price() : -std::numeric_limits<double>::infinity())
         < (b->hasMarket() ? b->market()->price() : -std::numeric_limits<double>::infinity());
}
bool BookStore::greater_price(const eris::SharedMember<Book> &a, const eris::SharedMember<Book> &b) {
    return (a->hasMarket() ? a->market()->price() : -std::numeric_limits<double>::infinity())
         > (b->hasMarket() ? b->market()->price() : -std::numeric_limits<double>::infinity());
}
bool BookStore::less_copies(const eris::SharedMember<Book> &a, const eris::SharedMember<Book> &b) {
    if (copies_cache_.count(a) != 0) copies_cache_[a] = a->copies();
    if (copies_cache_.count(b) != 0) copies_cache_[b] = b->copies();
    return copies_cache_[a] < copies_cache_[b];
}
bool BookStore::greater_copies(const eris::SharedMember<Book> &a, const eris::SharedMember<Book> &b) {
    if (copies_cache_.count(a) != 0) copies_cache_[a] = a->copies();
    if (copies_cache_.count(b) != 0) copies_cache_[b] = b->copies();
    return copies_cache_[a] > copies_cache_[b];
}

void BookStore::appendColumnsTo(Gtk::TreeView &v) const {
    appendCol(v, "ID", columns.id, 100);
    appendCol(v, "Position", columns.posstr, 150);
    if (not author_.ptr()) appendCol(v, "Author", columns.author, 100);
    appendCol(v, "Age", columns.age, 100);
    appendCol(v, "Quality", columns.quality, 100);
    appendCol(v, "Mkt?", columns.hasMarket, 50);
    appendCol(v, "Price", columns.price, 100);
    appendCol(v, "Rev.", columns.revenueLifetime, 100);
    appendCol(v, "Sales", columns.salesLifetime, 100);
    appendCol(v, "Copies", columns.copies, 100);
    appendCol(v, "Life", columns.lifetime, 100);
}

}}
