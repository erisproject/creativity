#include "creativity/gui/InfoWindow.hpp"
#include "creativity/gui/GUI.hpp"
#include <iomanip>
#include <gtkmm/box.h>
#include <gtkmm/treemodelcolumn.h>

using namespace eris;
using namespace creativity::state;

namespace creativity { namespace gui {

#define NEW_GRID(GRID) \
    grids_.emplace_back(); \
    Gtk::Grid &GRID = grids_.back(); \
    GRID.set_row_spacing(5); \
    GRID.set_column_spacing(5); \
    GRID.set_row_homogeneous(false); \
    GRID.set_column_homogeneous(true); \
    GRID.set_margin_top(10); \
    GRID.set_margin_left(10); \
    GRID.set_margin_bottom(10); \
    GRID.set_margin_right(10); \
    gpos = 0;

#define NEW_TAB_GRID(GRID, NB, TABNAME) \
    NEW_GRID(GRID); \
    NB.append_page(GRID, TABNAME);


// Creates a "value: data" row with labels stored in fields_ that can be updated later.  Unlike
// DATA_ROW, this takes a ROW and COL, and does not use gpos; the two labels are at (ROW,COL)
// and (ROW,COL+1)
#define DATA_AT(GRID, CODE, LABEL, ROW, COL) do { \
    fields_[CODE].first.set_markup(LABEL ":"); \
    fields_[CODE].first.set_alignment(1); \
    fields_[CODE].second.set_alignment(0); \
    GRID.attach(fields_[CODE].first, COL, ROW, 1, 1); \
    GRID.attach(fields_[CODE].second, COL+1, ROW, 1, 1); \
    } while (0)

/// Creates a "value: data" row with labels stored in fields_ that can be updated later
#define DATA_ROW(GRID, CODE, LABEL) do { \
    DATA_AT(GRID, CODE, LABEL, gpos, 0); \
    gpos++; \
    } while (0)

// Creates a single label at the given location, taking up the given number of rows/columns
#define LABEL_AT(GRID, LABEL, ROW, COL, WIDTH, HEIGHT, ALIGNMENT) do { \
    labels_.emplace_back(); Gtk::Label &a = labels_.back(); \
    a.set_alignment(ALIGNMENT); a.set_markup(LABEL); \
    GRID.attach(a, COL, ROW, WIDTH, HEIGHT); \
    } while (0)

// Creates a fixed "value: value" label pair on the given row in columns COL and COL+1.
#define LABELS_AT(GRID, LABEL, VALUE, ROW, COL) do { \
    LABEL_AT(GRID, LABEL ":", ROW, COL, 1, 1, 1); \
    LABEL_AT(GRID, VALUE,   ROW, COL+1, 1, 1, 0); \
    } while (0)

/// Creates a fixed "value: value" row
#define LABELS_ROW(GRID, LABEL, VALUE) do { \
    LABELS_AT(GRID, LABEL, VALUE, gpos, 0); \
    gpos++; \
    } while (0)

/// Creates a "value" row (which spans both columns)
#define COMMENT_ROW(GRID, COMMENT) do { \
    labels_.emplace_back(); \
    Gtk::Label &l = labels_.back(); \
    l.set_markup(COMMENT); \
    l.set_line_wrap(true); \
    l.set_max_width_chars(80); \
    l.set_margin_top(10); \
    l.set_margin_bottom(10); \
    GRID.attach(l, 0, gpos++, 2, 1); \
    } while (0)

#define BETA "<u><i>&#x3b2;</i></u>"
#define TIMES "&#xd7;"
#define PI "&#x3c0;"
#define GTREQ "&#x2265;"

InfoWindow::InfoWindow(std::shared_ptr<const State> state, std::shared_ptr<Gtk::Window> main_window,
        eris_id_t reader_id, std::function<void(eris::eris_id_t)> open_info_dialog)
    : reader{reader_id}, book{0}, open_info_dialog_{std::move(open_info_dialog)}
{
    if (state->readers.count(reader) == 0)
        throw std::out_of_range("InfoWindow() called with invalid reader id");

    initWindow(*main_window);

    set_title("Reader/author details (" + std::to_string(reader) + ")");

    nbs_.emplace_back();
    Gtk::Notebook &tabs = nbs_.back();

    int gpos;

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
    LABELS_ROW(grid_quality, "Dependent variable", "<i>quality</i>");
    DATA_ROW(grid_quality, "q_n", "n");
    std::vector<std::string> q_vars{{"constant", "I(firstBook)", "prevBooks", "age", "price", "price" TIMES "age", "copiesSold"}};
    for (size_t i = 0; i < q_vars.size(); i++)
        DATA_ROW(grid_quality, "q_" + std::to_string(i), BETA "[" + q_vars[i] + "]");

    NEW_TAB_GRID(grid_profit, beliefs, "Profit");
    LABELS_ROW(grid_profit, "Dependent variable", "<i>lifetimeProfit</i>");
    DATA_ROW(grid_profit, "p_n", "n");
    std::vector<std::string> p_vars{{"constant", "quality<sup>D</sup>", "I(firstBook)", "previousBooks", "marketBooks"}};
    for (size_t i = 0; i < p_vars.size(); i++)
        DATA_ROW(grid_profit, "p_" + std::to_string(i), BETA "[" + p_vars[i] + "]");
    DATA_ROW(grid_profit, "_p_draws", "# successful draws");
    DATA_ROW(grid_profit, "_p_discards", "# discarded draws");

    NEW_TAB_GRID(grid_demand, beliefs, "Demand");
    COMMENT_ROW(grid_demand, "Note: this regression is for single-period demand");
    LABELS_ROW(grid_demand, "Dependent variable", "<i>quantityDemanded</i>");
    DATA_ROW(grid_demand, "d_n", "n");
    std::vector<std::string> d_vars{{"constant", "price<sup>D</sup>", "quality<sup>D</sup>", "prevSales", "age", "I(onlyBook)", "otherBooks", "marketBooks"}};
    for (size_t i = 0; i < d_vars.size(); i++)
        DATA_ROW(grid_demand, "d_" + std::to_string(i), BETA "[" + d_vars[i] + "]");
    DATA_ROW(grid_demand, "_d_draws", "# successful draws");
    DATA_ROW(grid_demand, "_d_discards", "# discarded draws");

    NEW_TAB_GRID(grid_pstream, beliefs, "Profit Stream");

    LABEL_AT(grid_pstream, "Dependent variable: <i>profitRemaining</i>", 0, 0, 2*Reader::profit_stream_ages.size(), 1, 0.5);

    int col = 0;
    for (unsigned long a : Reader::profit_stream_ages) {
        std::string head_label = "<u>age " GTREQ " " + std::to_string(a) + "</u>";
        std::string code_prefix = "ps" + std::to_string(a) + "_";
        LABEL_AT(grid_pstream, head_label, 1, col, 2, 1, 0.5);

        DATA_AT(grid_pstream, code_prefix + "n", "n", 2, col);
        for (unsigned long j = 0; j < a; j++) {
            DATA_AT(grid_pstream, code_prefix + std::to_string(j), BETA "[" PI "<sub>" + std::to_string(j) + "</sub>]", 3+j, col);
        }
        col += 2;
    }

    NEW_TAB_GRID(grid_pextrap, beliefs, "Profit (extrap.)");
    LABELS_ROW(grid_pextrap, "Dependent variable", "<i>lifetimeProfit</i>");
    DATA_ROW(grid_pextrap, "pe_n", "n");
    for (size_t i = 0; i < p_vars.size(); i++)
        DATA_ROW(grid_pextrap, "pe_" + std::to_string(i), BETA "[" + p_vars[i] + "]");
    COMMENT_ROW(grid_pextrap, "<i>NB: This is the same model as the Profit belief, but its data also includes extrapolated values for "
            "still-on-market books using ProfitStream beliefs, while Profit beliefs only include books once they leave the market.</i>");

    bk_model_ = BookStore::create(state, reader);
    bk_tree_.set_model(bk_model_);
    bk_model_->appendColumnsTo(bk_tree_);
    bk_tree_.set_fixed_height_mode(true);
    bk_model_->set_sort_column(bk_model_->columns.id, Gtk::SortType::SORT_DESCENDING);
    bk_tree_.signal_row_activated().connect([this] (const Gtk::TreeModel::Path &path, Gtk::TreeViewColumn*) -> void {
        open_info_dialog_(bk_model_->member(path).id);
    });
    swins_.emplace_back();
    swins_.back().add(bk_tree_);
    swins_.back().set_policy(Gtk::POLICY_NEVER, Gtk::POLICY_AUTOMATIC);
    tabs.append_page(swins_.back(), "Authored Books");

    add(tabs);

    refresh(state);

    show_all();
}

InfoWindow::InfoWindow(std::shared_ptr<const State> state, std::shared_ptr<Gtk::Window> main_window, unsigned long book_id)
    : reader{0}, book{book_id}
{
    initWindow(*main_window);

    set_title("Book details (" + std::to_string(book) + ")");

    int gpos;

    NEW_GRID(grid_status);

    DATA_ROW(grid_status, "id", "ID");
    DATA_ROW(grid_status, "position", "Position");
    DATA_ROW(grid_status, "market", "Market");
    DATA_ROW(grid_status, "price", "Price");
    DATA_ROW(grid_status, "quality", "Quality");
    DATA_ROW(grid_status, "age", "Age");
    DATA_ROW(grid_status, "revenue", "Revenue (lifetime)");
    DATA_ROW(grid_status, "revenueLast", "Revenue (current)");
    DATA_ROW(grid_status, "sales", "Copies sold (lifetime)");
    DATA_ROW(grid_status, "salesLast", "Copies sold (current)");
    DATA_ROW(grid_status, "copies", "Copies in world");
    DATA_ROW(grid_status, "author", "Author ID");
    override_background_color(Gdk::RGBA{"white"});

    add(grid_status);

    refresh(state);

    show_all();
}

void InfoWindow::initWindow(Gtk::Window &parent) {
    set_transient_for(parent);
//    set_type_hint(Gdk::WINDOW_TYPE_HINT_UTILITY);

    set_resizable(true);
    property_destroy_with_parent() = true;
    set_deletable(true);
}

void InfoWindow::refresh(std::shared_ptr<const State> state) {
    if (state->t == t_) return;

    ERIS_DBG("refreshing");

    if (reader) {
        auto &r = state->readers.at(reader);
        updateValue("id", r.id);
        updateValue("position", GUI::pos_to_string(r.position));
        updateValue("utility", r.u);
        updateValue("uLife", r.u_lifetime);
        updateValue("books", r.library.size());
        updateValue("booksNew", r.new_books.size());
        updateValue("booksWritten", r.wrote.size());
        updateValue("bookLast", r.wrote.empty()
                ? "(never written)"
                : std::to_string(state->books.at(r.wrote.back()).age));

#define UPDATE_LIN(PREFIX, VAR) \
        updateValue(PREFIX + std::string("n"), VAR.n()); \
        for (size_t i = 0; i < VAR.K(); i++) \
            updateValue(PREFIX + std::to_string(i), VAR.beta()[i]);
#define UPDATE_LIN_RB(PREFIX, BELIEF) UPDATE_LIN(PREFIX, r.BELIEF)

        UPDATE_LIN_RB("q_", quality);
        UPDATE_LIN_RB("p_", profit);
        UPDATE_LIN_RB("d_", demand);
        UPDATE_LIN_RB("pe_", profit_extrap);

        updateValue("_p_draws", r.profit.draw_success_cumulative);
        updateValue("_p_discards", r.profit.draw_discards_cumulative);
        updateValue("_d_draws", r.demand.draw_success_cumulative);
        updateValue("_d_discards", r.demand.draw_discards_cumulative);

        for (unsigned long a : Reader::profit_stream_ages) {
            std::string code_prefix = "ps" + std::to_string(a) + "_";

            if (r.profit_stream.count(a)) {
                auto &psi = r.profit_stream.at(a);
                UPDATE_LIN(code_prefix, psi);
            }
            else {
                updateValue(code_prefix + "n", "0");
                for (size_t j = 0; j < a; j++)
                    updateValue(code_prefix + std::to_string(j), "");
            }
        }

        // Update the books tree
        if (not initial_refresh_) {
            // Preserve the current sort column/order, if any
            int sort_col;
            Gtk::SortType sort_order;
            bool resort = bk_model_->get_sort_column_id(sort_col, sort_order);
            bk_model_ = BookStore::create(state, reader);
            if (resort) bk_model_->set_sort_column(sort_col, sort_order);
            bk_tree_.set_model(bk_model_);
        }
    }
    else {
        if (state->books.count(book)) {
            set_title("Book details (" + std::to_string(book) + ")");
            auto &b = state->books.at(book);
            updateValue("id", b.id);
            updateValue("author", b.author);
            updateValue("position", GUI::pos_to_string(b.position));
            updateValue("market", b.market ? "yes" : "no");
            updateValue("price", b.price);
            updateValue("quality", b.quality);
            updateValue("age", b.age);
            updateValue("revenue", b.revenue_lifetime);
            updateValue("revenueLast", b.revenue);
            updateValue("sales", b.sales_lifetime);
            updateValue("salesLast", b.sales);
            updateValue("copies", b.copies);
        }
        else {
            // If the user navigates back in time to a period where the book doesn't exist, set
            // everything to N/A and change the title.
            set_title("Book details (" + std::to_string(book) + "): not yet created!");
            for (auto &key : {"id", "author", "position", "market", "price", "quality", "age", "revenue", "revenueLast", "sales", "salesLast", "copies"}) {
                updateValue(key, "N/A");
            }
        }
    }

    if (initial_refresh_) initial_refresh_ = false;
}

void InfoWindow::updateValue(const std::string &code, const std::string &val) {
    fields_[code].second.set_text(val);
}
void InfoWindow::updateValue(const std::string &code, unsigned long val) {
    updateValue(code, std::to_string(val));
}
void InfoWindow::updateValue(const std::string &code, long val) {
    updateValue(code, std::to_string(val));
}
void InfoWindow::updateValue(const std::string &code, double val) {
    updateValue(code, std::to_string(val));
}



} }
