#include "creativity/gui/LibraryStore.hpp"
#include "creativity/state/BookState.hpp"

namespace creativity { namespace gui {

using namespace eris;
using namespace creativity::state;
using namespace std::placeholders;

LibraryStore::LibraryStore(std::shared_ptr<const State> &&state, eris_id_t reader)
    : Glib::ObjectBase(typeid(LibraryStore)),
    BookStore(std::move(state), std::unique_ptr<BookStore::ColRec>(new LibraryStore::ColRec())),
    reader_(std::move(reader))
{
    // Add the reader's library:
    for (auto &l : state_->readers.at(reader_).library) {
        if (not l.second.wrote()) members_.emplace_back(state_->books.at(l.first));
    }
}

Glib::RefPtr<LibraryStore> LibraryStore::create(std::shared_ptr<const State> state, eris_id_t reader) {
    return Glib::RefPtr<LibraryStore>(new LibraryStore(std::move(state), reader));
}

void LibraryStore::get_value_vfunc(const iterator &iter, int column, Glib::ValueBase &value) const {
    if (iter.get_stamp() != stamp_ or column > get_n_columns_vfunc()) return;

    const LibraryStore::ColRec &cols = static_cast<const LibraryStore::ColRec&>(*columns);

// If statement for field with accessor A into the reader state library element:
#define IFCOL(F, A) if (column == cols.F.index()) copy_value_(value, \
        state_->readers.at(reader_).library.at( \
                members_.at((size_t) iter.gobj()->user_data).get().id \
            ).A)
    IFCOL(reader_quality, quality);
    else IFCOL(reader_pirated, pirated());
    else IFCOL(reader_purchased, purchased());
    else IFCOL(reader_acquired, acquired);
#undef IFCOL
    else {
        // Otherwise it's a book property, pass through
        BookStore::get_value_vfunc(iter, column, value);
    }
}

void LibraryStore::set_sort_column_id_vfunc(int sort_column_id, Gtk::SortType order) {
    bool ascending = (order == Gtk::SORT_ASCENDING);
    std::function<bool(const BookState &a, const BookState &b)> compare;

    const LibraryStore::ColRec &cols = static_cast<const LibraryStore::ColRec&>(*columns);
#define IF_COL(COL) if (sort_column_id == cols.COL.index()) \
        compare = std::bind(ascending ? &LibraryStore::less_##COL : &LibraryStore::greater_##COL, this, _1, _2)
    IF_COL(reader_quality);
    else IF_COL(reader_pirated);
    else IF_COL(reader_purchased);
    else IF_COL(reader_acquired);
#undef IF_COL
    else return BookStore::set_sort_column_id_vfunc(sort_column_id, order);

    sort_members(compare, sort_column_id, order);
}

#define LESS_GREATER_METHODS(COL, MEMBER) \
bool LibraryStore::less_##COL(const state::BookState &a, const state::BookState &b) { \
    auto &lib = state_->readers.at(reader_).library; \
    return lib.at(a.id).MEMBER < lib.at(b.id).MEMBER; \
} \
bool LibraryStore::greater_##COL(const state::BookState &a, const state::BookState &b) { \
    auto &lib = state_->readers.at(reader_).library; \
    return lib.at(a.id).MEMBER > lib.at(b.id).MEMBER; \
}
LESS_GREATER_METHODS(reader_quality, quality)
LESS_GREATER_METHODS(reader_purchased, purchased())
LESS_GREATER_METHODS(reader_pirated, pirated())
LESS_GREATER_METHODS(reader_acquired, acquired)
#undef LESS_GREATER_METHODS

void LibraryStore::appendColumnsTo(Gtk::TreeView &v) const {
    auto &col = static_cast<LibraryStore::ColRec&>(*columns);
    appendCol(v, "ID", col.id, 50);
    appendCol(v, "Position", col.pos_str, 110);
    appendCol(v, "Author", col.author, 80);
    appendCol(v, "Created", col.created, 85);
    appendCol(v, "Acquired", col.reader_acquired, 85);
    appendCol(v, "Bought", col.reader_purchased, 65);
    appendCol(v, "Pirated", col.reader_pirated, 65);
    appendCol(v, "Quality", col.reader_quality, 85);
    appendCol(v, "Base Q.", col.quality, 85);
}

}}
