#include "creativity/GUIInfoWindow.hpp"
#include "creativity/Book.hpp"
#include "creativity/Reader.hpp"
#include <iomanip>
#include <gtkmm/box.h>

namespace creativity {

using namespace eris;

#define NEW_GRID(GRID) \
    grids_.emplace_back(); \
    Gtk::Grid &GRID = grids_.back(); \
    GRID.set_row_spacing(5); \
    GRID.set_column_spacing(5); \
    GRID.set_row_homogeneous(false); \
    GRID.set_column_homogeneous(true); \
    gpos = 0;

#define NEW_TAB_GRID(GRID, NB, tabname) \
    NEW_GRID(GRID); \
    NB.append_page(GRID, tabname);

/// Creates a "value: data" row with labels stored in fields_ that can be updated later
#define DATA_ROW(grid, code, label) do { \
    fields_[code].first.set_markup(label ":"); \
    fields_[code].first.set_alignment(1); \
    fields_[code].second.set_alignment(0); \
    grid.attach(fields_[code].first, 0, gpos, 1, 1); \
    grid.attach(fields_[code].second, 1, gpos++, 1, 1); \
    } while (0)

/// Creates a fixed "value: value" row
#define LABEL_ROW(grid, label, value) do { \
    labels_.emplace_back(); Gtk::Label &a = labels_.back(); \
    labels_.emplace_back(); Gtk::Label &b = labels_.back(); \
    a.set_alignment(1); a.set_markup(label ":"); \
    b.set_alignment(0); b.set_markup(value); \
    grid.attach(a, 0, gpos, 1, 1); \
    grid.attach(b, 1, gpos++, 1, 1); \
    } while (0)


/// Creates a "value" row (which spans both columns)
#define COMMENT_ROW(grid, comment) do { \
    labels_.emplace_back(); \
    Gtk::Label &l = labels_.back(); \
    l.set_markup(comment); \
    grid.attach(l, 0, gpos++, 2, 1); \
    } while (0)

#define BETA "<u><i>&#x3b2;</i></u>"
#define TIMES "&#xd7;"
#define PI "&#x3c0;"


GUIInfoWindow::GUIInfoWindow(eris::SharedMember<Reader> rdr, std::shared_ptr<Gtk::Window> main_window)
    : reader{rdr}
{
    set_transient_for(*main_window);
//    set_type_hint(Gdk::WINDOW_TYPE_HINT_UTILITY);
    set_title("Reader/author details");
    set_resizable(true);
    property_destroy_with_parent() = true;
    set_deletable(true);

    int gpos;

    nbs_.emplace_back();
    Gtk::Notebook &tabs = nbs_.back();

    NEW_TAB_GRID(grid_status, tabs, "Status");
    DATA_ROW(grid_status, "id", "ID");
    DATA_ROW(grid_status, "position", "Position");
    DATA_ROW(grid_status, "utility", "Utility");
    DATA_ROW(grid_status, "uLife", "Lifetime utility");
    DATA_ROW(grid_status, "books", "Books owned");
    DATA_ROW(grid_status, "booksNew", "New books");
    DATA_ROW(grid_status, "booksWritten", "Books written");
    DATA_ROW(grid_status, "bookLast", "Latest book age");

    nbs_.emplace_back();
    Gtk::Notebook &beliefs = nbs_.back();
    tabs.append_page(beliefs, "Beliefs");

    NEW_TAB_GRID(grid_quality, beliefs, "Quality");
    LABEL_ROW(grid_quality, "Dependent variable", "<i>quality</i>");
    std::vector<std::string> q_vars{{"constant", "I(firstBook)", "prevBooks", "age", "price", "price" TIMES "age", "copiesSold"}};
    for (int i = 0; i < reader->qualityBelief().K(); i++)
        DATA_ROW(grid_quality, "q_" + std::to_string(i), BETA "[" + q_vars[i] + "]");

    NEW_TAB_GRID(grid_profit, beliefs, "Profit");
    LABEL_ROW(grid_profit, "Dependent variable", "<i>lifetimeProfit</i>");
    std::vector<std::string> p_vars{{"constant", "quality<sup>D</sup>", "I(firstBook)", "previousBooks", "marketBooks"}};
    for (int i = 0; i < reader->profitBelief().K(); i++)
        DATA_ROW(grid_profit, "p_" + std::to_string(i), BETA "[" + p_vars[i] + "]");

    NEW_TAB_GRID(grid_demand, beliefs, "Demand");
    LABEL_ROW(grid_demand, "Dependent variable", "<i>quantityDemanded</i>");
    COMMENT_ROW(grid_demand, "Note: this regression is for single-period demand");
    std::vector<std::string> d_vars{{"constant", "price<sup>D</sup>", "quality<sup>D</sup>", "prevSales", "age", "I(onlyBook)", "otherBooks", "marketBooks"}};
    for (int i = 0; i < reader->demandBelief().K(); i++)
        DATA_ROW(grid_demand, "d_" + std::to_string(i), BETA "[" + d_vars[i] + "]");

    NEW_TAB_GRID(grid_pstream, beliefs, "Profit stream");
    LABEL_ROW(grid_pstream, "Dependent variable", "<i>profitRemaining</i>");
    for (int i = 0; i < reader->profitStreamBelief().K(); i++) {
        DATA_ROW(grid_pstream, "ps_" + std::to_string(i),
                BETA "[" + (i == 0 ? "constant" : "I(age=" + std::to_string(i) + ") " TIMES " " PI "<sub>" + std::to_string(i-1) + "</sub>") + "]");
    }

    NEW_TAB_GRID(grid_pextrap, beliefs, "Profit (extrap.)");
    LABEL_ROW(grid_pextrap, "Dependent variable", "<i>lifetimeProfit</i>");
    for (int i = 0; i < reader->profitExtrapBelief().K(); i++)
        DATA_ROW(grid_pextrap, "pe_" + std::to_string(i), BETA "[" + p_vars[i] + "]");
    COMMENT_ROW(grid_pextrap, "<i>NB: regression includes still-on-market books using profit stream prediction</i>");

    gpos = 0;

    // Call refresh to actually set all the above values to the right thing:
    refresh();

    add(tabs);

    show_all();
}

GUIInfoWindow::GUIInfoWindow(eris::SharedMember<Book> bk, std::shared_ptr<Gtk::Window> main_window)
    : book{bk}
{
    set_transient_for(*main_window);
//    set_type_hint(Gdk::WINDOW_TYPE_HINT_UTILITY);
    set_title("Book details");
    set_resizable(true);
    property_destroy_with_parent() = true;
    set_deletable(true);

    int gpos;

    NEW_GRID(grid_status);

    DATA_ROW(grid_status, "id", "ID");
    DATA_ROW(grid_status, "position", "Position");
    DATA_ROW(grid_status, "market", "Market ID");
    DATA_ROW(grid_status, "price", "Price");
    DATA_ROW(grid_status, "quality", "Quality");
    DATA_ROW(grid_status, "age", "Age");
    DATA_ROW(grid_status, "revenue", "Revenue (lifetime)");
    DATA_ROW(grid_status, "revenueLast", "Revenue (last)");
    DATA_ROW(grid_status, "sales", "Copies sold (lifetime)");
    DATA_ROW(grid_status, "salesLast", "Copies sold (last)");
    DATA_ROW(grid_status, "copies", "Copies in world");
    DATA_ROW(grid_status, "author", "Author ID");

    refresh();

    override_background_color(Gdk::RGBA{"white"});
    add(grid_status);

    show_all();
}

void GUIInfoWindow::refresh() {
    if (reader) {
        updateValue("id", reader->id());
        updateValue("position", pos_to_string(reader->position()));
        updateValue("utility", reader->u());
        updateValue("uLife", reader->uLifetime());
        updateValue("books", reader->library().size());
        updateValue("booksNew", reader->newBooks().size());
        updateValue("booksWritten", reader->wrote().size());
        updateValue("bookLast", reader->wrote().empty()
                ? "(never written)"
                : std::to_string(reader->wrote().back()->age()));

#define UPDATE_LIN(PREFIX, BELIEF) \
        for (int i = 0; i < reader->BELIEF().K(); i++) \
            updateValue(PREFIX + std::to_string(i), reader->BELIEF().betaPrior()[i]);

        UPDATE_LIN("q_", qualityBelief);
        UPDATE_LIN("p_", profitBelief);
        UPDATE_LIN("d_", demandBelief);
        UPDATE_LIN("ps_", profitStreamBelief);
        UPDATE_LIN("pe_", profitExtrapBelief);
    }
    else {
        updateValue("id", book->id());
        updateValue("author", book->author()->id());
        updateValue("position", pos_to_string(book->position()));
        if (book->hasMarket()) {
            updateValue("market", book->market()->id());
            updateValue("price", book->market()->price());
        }
        else {
            updateValue("market", "(not on market)");
            updateValue("price", "(not on market)");
        }
        updateValue("quality", book->quality());
        updateValue("age", book->age());
        updateValue("revenue", book->lifeRevenue());
        updateValue("revenueLast", book->currRevenue());
        updateValue("sales", book->lifeSales());
        updateValue("salesLast", book->currSales());

        unsigned long copies = 0;
        for (auto &r : book->simulation()->agents<Reader>()) {
            if (r->library().count(book))
                copies++;
        }
        updateValue("copies", copies);
    }
}

void GUIInfoWindow::updateValue(const std::string &code, const std::string &val) {
    fields_[code].second.set_text(val);
}
void GUIInfoWindow::updateValue(const std::string &code, unsigned long val) {
    updateValue(code, std::to_string(val));
}
void GUIInfoWindow::updateValue(const std::string &code, long val) {
    updateValue(code, std::to_string(val));
}
void GUIInfoWindow::updateValue(const std::string &code, double val) {
    updateValue(code, std::to_string(val));
}

std::string GUIInfoWindow::pos_to_string(const eris::Position &pos) {
    std::ostringstream opos;
    opos << "(" << std::setw(7) << std::showpoint;
    for (size_t i = 0; i < pos.dimensions; i++) {
        if (i > 0) opos << ",";
        opos << pos[i];
    }
    opos << ")";
    return opos.str();
}

}
