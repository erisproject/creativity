#include "creativity/gui/BookStore.hpp"
#include "creativity/gui/GUI.hpp"
#include <gtkmm/treepath.h>
#include <algorithm>

using namespace eris;
using namespace std::placeholders;

namespace creativity { namespace gui {
    
BookStore::BookStore(std::shared_ptr<Simulation> &&sim, SharedMember<Reader> &&author)
    : Glib::ObjectBase(typeid(BookStore)), MemberStore(std::move(sim)),
    author_{std::move(author)}
{}

Glib::RefPtr<BookStore> BookStore::create(std::shared_ptr<Simulation> sim, SharedMember<Reader> author) {
    return Glib::RefPtr<BookStore>(new BookStore(std::move(sim), std::move(author)));
}

std::vector<SharedMember<Book>> BookStore::resync_add() {
    copies_cache_.clear();

    std::vector<SharedMember<Book>> new_books;
    // NB: check .ptr rather than ->id() because the author might have been removed from the simulation,
    // but we are still responsible for that author's books and don't want to revert to all simulation books.
    if (author_.ptr()) {
        for (auto &b : author_->wrote()) {
            if (b > max_id_) new_books.push_back(b);
        }

        return new_books;
    }

    // else global book list:
    return sim_->goods<Book>([this](const Book &b) -> bool { return b.id() > max_id_; });
}

int BookStore::get_n_columns_vfunc() const {
    return columns.size();
}

GType BookStore::get_column_type_vfunc(int index) const {
    if ((unsigned int) index >= columns.size()) throw std::out_of_range("Invalid column index accessed");
    return columns.types()[index];
}

void BookStore::get_value_vfunc(const iterator &iter, int column, Glib::ValueBase &value) const {
    if (iter.get_stamp() != stamp_ or column > get_n_columns_vfunc()) return;

    auto &b = members_.at((size_t) iter.gobj()->user_data);
    if (column == columns.id.index() || column == columns.author.index() || column == columns.market.index()) {
        Glib::Value<eris_id_t> v;
        v.init(v.value_type());
        v.set(  column == columns.id.index() ? b->id() :
                column == columns.author.index() ? b->author()->id() :
                b->hasMarket() ? b->market()->id() : 0
             );
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
        Glib::Value<bool> v;
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
        value.init(v.gobj());
    }
    else if (column == columns.age.index() or column == columns.sales.index() or column == columns.salesLifetime.index()
            or column == columns.lifetime.index()) {
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

void BookStore::set_sort_column_id_vfunc(int sort_column_id, Gtk::SortType order) {
    bool ascending = (order == Gtk::SORT_ASCENDING);
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

    sort_members(compare, sort_column_id, order);
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
    appendCol(v, "ID", columns.id, 50);
    appendCol(v, "Position", columns.posstr, 100);
    if (not author_.ptr()) appendCol(v, "Author", columns.author, 80);
    appendCol(v, "Age", columns.age, 65);
    appendCol(v, "Quality", columns.quality, 85);
    appendCol(v, "Mkt?", columns.hasMarket, 65);
    appendCol(v, "Price", columns.price, 75);
    appendCol(v, "Rev.", columns.revenueLifetime, 75);
    appendCol(v, "Sales", columns.salesLifetime, 75);
    appendCol(v, "Copies", columns.copies, 80);
    appendCol(v, "Life", columns.lifetime, 65);
}

}}
