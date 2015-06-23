#include "creativity/gui/InfoWindow.hpp"
#include "creativity/gui/GUI.hpp"
#include <iomanip>
#include <gtkmm/box.h>
#include <gtkmm/treemodelcolumn.h>
#include <Eigen/QR>

using namespace eris;
using namespace creativity::state;
using namespace eris::belief;
using namespace Eigen;

namespace creativity { namespace gui {

Gtk::Grid& InfoWindow::new_grid() {
    grids_.emplace_back();
    Gtk::Grid &g = grids_.back();
    g.set_row_spacing(5);
    g.set_column_spacing(5);
    g.set_row_homogeneous(false);
    g.set_column_homogeneous(false);
    g.set_margin_top(10);
    g.set_margin_left(10);
    g.set_margin_bottom(10);
    g.set_margin_right(10);
    gpos_ = 0;
    return g;
}

Gtk::Grid& InfoWindow::new_tab_grid(Gtk::Notebook &notebook, const std::string &title) {
    auto &g = new_grid();
    notebook.append_page(g, title);
    return g;
}

void InfoWindow::data_at(Gtk::Grid &grid, const std::string &code, const std::string &value_name, int row, int col,
        double label_align, double value_align) {
    fields_[code].first.set_markup(value_name + ":");
    fields_[code].first.set_alignment(label_align);
    fields_[code].second.set_alignment(value_align);
    grid.attach(fields_[code].first, col, row, 1, 1);
    grid.attach(fields_[code].second, col+1, row, 1, 1);
}

void InfoWindow::data_append(Gtk::Grid &grid, const std::string &code, const std::string &value_name,
        double label_align, double value_align) {
    data_at(grid, code, value_name, gpos_, 0, label_align, value_align);
    gpos_++;
}

void InfoWindow::matrix_at(Gtk::Grid &grid, const std::string &code, const std::string &name, int row, int col, int nrows, int ncols) {
    // Matrix label:
    label_at(grid, name, 1.0, row, col, 1, 1);

    // Column labels:
    for (int j = 1; j <= ncols; j++) label_at(grid, "[," + std::to_string(j) + "]", 1.0, row, col+j, 1, 1);

    auto &labels = matrix_[code];
    for (int i = 1; i <= nrows; i++) {
        // Row label:
        label_at(grid, "[" + std::to_string(i) + ",]", 1.0, row+i, col, 1, 1);

        // Matrix element labels:
        for (int j = 1; j <= ncols; j++) {
            labels.emplace_back(new Gtk::Label);
            Gtk::Label &l = *labels.back();
            l.set_alignment(1);
            grid.attach(l, col+j, row+i, 1, 1);
        }
    }
}


void InfoWindow::label_at(Gtk::Grid &grid, const std::string &label, double alignment, int row, int col, int width, int height) {
    labels_.emplace_back();
    Gtk::Label &a = labels_.back();
    a.set_alignment(alignment);
    a.set_markup(label);
    grid.attach(a, col, row, width, height);
}

void InfoWindow::labels_at(Gtk::Grid &grid, const std::string &label1, const std::string &label2, int row, int col) {
    label_at(grid, label1 + ":", 1.0, row, col, 1, 1);
    label_at(grid, label2, 0.0, row, col+1, 1, 1);
}
void InfoWindow::labels_append(Gtk::Grid &grid, const std::string &label1, const std::string &label2) {
    labels_at(grid, label1, label2, gpos_++, 0);
}
void InfoWindow::comment_append(Gtk::Grid &grid, const std::string &comment, int cols) {
    labels_.emplace_back();
    auto &c = labels_.back();
    c.set_markup(comment);
    c.set_line_wrap(true);
    c.set_max_width_chars(80);
    c.set_margin_top(10);
    c.set_margin_bottom(10);
    grid.attach(c, 0, gpos_++, cols, 1);
}


#define BETA "<i>β</i>"

InfoWindow::InfoWindow(std::shared_ptr<const State> state, std::shared_ptr<Gtk::Window> main_window,
        eris_id_t reader_id, std::function<void(eris::eris_id_t)> open_info_dialog)
    : reader{reader_id}, book{0}, open_info_dialog_{std::move(open_info_dialog)}
{
    if (state->readers.count(reader) == 0)
        throw std::out_of_range("InfoWindow() called with invalid reader id");

    auto rdr = state->readers.at(reader);

    initWindow(*main_window);

    set_title("Reader/author details (" + std::to_string(reader) + ")");

    nbs_.emplace_back();
    Gtk::Notebook &tabs = nbs_.back();

    auto &grid_status = new_tab_grid(tabs, "Status");
    data_append(grid_status, "id", "ID");
    data_append(grid_status, "position", "Position");
    data_append(grid_status, "utility", "Utility");
    data_append(grid_status, "uLife", "Lifetime utility");
    data_append(grid_status, "books", "Books owned");
    data_append(grid_status, "booksPurchased", "Books purchased");
    data_append(grid_status, "booksPirated", "Books pirated");
    data_append(grid_status, "booksNew", "New books");
    data_append(grid_status, "booksNewPurchased", "New books (purchased)");
    data_append(grid_status, "booksNewPirated", "New books (pirated)");
    data_append(grid_status, "booksWritten", "Books written");
    data_append(grid_status, "bookLast", "Latest book age");
    data_append(grid_status, "numFriends", "# Friends");

    nbs_.emplace_back();
    Gtk::Notebook &beliefs = nbs_.back();
    tabs.append_page(beliefs, "Beliefs");

    auto &grid_quality = new_tab_grid(beliefs, "Quality");
    labels_append(grid_quality, "Dependent variable", "<i>quality</i>");
    data_append(grid_quality, "q_n", "n");
    data_append(grid_quality, "q_s2", "s<sup>2</sup>");
    std::vector<std::string> q_vars(rdr.quality.names());
    for (size_t i = 0; i < q_vars.size(); i++)
        data_append(grid_quality, "q_" + std::to_string(i), BETA "[" + q_vars[i] + "]", 1, 1);
    matrix_at(grid_quality, "q_V", "s<sup>2</sup><b>V</b>", 2, 2, q_vars.size(), q_vars.size());

    auto &grid_profit = new_tab_grid(beliefs, "Profit");
    labels_append(grid_profit, "Dependent variable", "<i>lifetimeProfit</i>");
    data_append(grid_profit, "p_n", "n");
    data_append(grid_profit, "p_s2", "s<sup>2</sup>");
    std::vector<std::string> p_vars(rdr.profit.names());
    for (size_t i = 0; i < p_vars.size(); i++)
        data_append(grid_profit, "p_" + std::to_string(i), BETA "[" + p_vars[i] + "]");
    data_append(grid_profit, "_p_draws", "# successful draws");
    data_append(grid_profit, "_p_discards", "# discarded draws");
    data_append(grid_profit, "_p_drawtype", "Draw method");
    matrix_at(grid_profit, "p_V", "s<sup>2</sup><b>V</b>", 2, 2, p_vars.size(), p_vars.size());

    auto &grid_demand = new_tab_grid(beliefs, "Demand");
    labels_append(grid_demand, "Dependent variable", "<i>quantityDemanded</i>");
    data_append(grid_demand, "d_n", "n");
    data_append(grid_demand, "d_s2", "s<sup>2</sup>");
    std::vector<std::string> d_vars(rdr.demand.names());
    for (size_t i = 0; i < d_vars.size(); i++)
        data_append(grid_demand, "d_" + std::to_string(i), BETA "[" + d_vars[i] + "]");
    data_append(grid_demand, "_d_draws", "# successful draws");
    data_append(grid_demand, "_d_discards", "# discarded draws");
    data_append(grid_demand, "_d_drawtype", "Draw method");
    matrix_at(grid_demand, "d_V", "s<sup>2</sup><b>V</b>", 2, 2, d_vars.size(), d_vars.size());
    comment_append(grid_demand, "<i>NB: This regression is for single-period demand.</i>", p_vars.size() + 3);

    nbs_.emplace_back();
    Gtk::Notebook &pstream_tabs = nbs_.back();
    beliefs.append_page(pstream_tabs, "Profit Stream");


    for (unsigned long a : Reader::profit_stream_ages) {
        auto &grid_pstream = new_tab_grid(pstream_tabs, "age ≥ " + std::to_string(a));
        labels_append(grid_pstream, "Dependent variable", "<i>profitRemaining</i>");
        std::string code_prefix = "ps" + std::to_string(a) + "_";

        data_append(grid_pstream, code_prefix + "n", "n");
        data_append(grid_pstream, code_prefix + "s2", "s<sup>2</sup>");
        for (unsigned long j = 0; j < a; j++) {
            data_append(grid_pstream, code_prefix + std::to_string(j), BETA "[π<sub>" + std::to_string(j) + "</sub>]");
        }
        matrix_at(grid_pstream, code_prefix + "V", "s<sup>2</sup><b>V</b>", 2, 2, a, a);
    }

    auto &grid_pextrap = new_tab_grid(beliefs, "Profit (extrap.)");
    labels_append(grid_pextrap, "Dependent variable", "<i>lifetimeProfit</i>");
    data_append(grid_pextrap, "pe_n", "n");
    data_append(grid_pextrap, "pe_s2", "s<sup>2</sup>");
    for (size_t i = 0; i < p_vars.size(); i++)
        data_append(grid_pextrap, "pe_" + std::to_string(i), BETA "[" + p_vars[i] + "]");
    data_append(grid_pextrap, "_pe_draws", "# successful draws");
    data_append(grid_pextrap, "_pe_discards", "# discarded draws");
    data_append(grid_pextrap, "_pe_drawtype", "Draw method");
    matrix_at(grid_pextrap, "pe_V", "s<sup>2</sup><b>V</b>", 2, 2, p_vars.size(), p_vars.size());
    comment_append(grid_pextrap, "<i>NB: This is the same model as the Profit belief, but its data also includes extrapolated values for "
            "still-on-market books using ProfitStream beliefs, while Profit beliefs only include books once they leave the market.</i>",
            p_vars.size() + 3);

    bk_authored_model_ = BookStore::create(state, reader);
    bk_authored_tree_.set_model(bk_authored_model_);
    bk_authored_model_->appendColumnsTo(bk_authored_tree_);
    bk_authored_tree_.set_fixed_height_mode(true);
    bk_authored_model_->set_sort_column(bk_authored_model_->columns->id, Gtk::SortType::SORT_DESCENDING);
    bk_authored_tree_.signal_row_activated().connect([this] (const Gtk::TreeModel::Path &path, Gtk::TreeViewColumn*) -> void {
        open_info_dialog_(bk_authored_model_->member(path).id);
    });
    swins_.emplace_back();
    swins_.back().add(bk_authored_tree_);
    swins_.back().set_policy(Gtk::POLICY_NEVER, Gtk::POLICY_AUTOMATIC);
    fields_["authored_books_tab"].second.set_text("Authored Books");
    tabs.append_page(swins_.back(), fields_["authored_books_tab"].second);

    bk_library_model_ = LibraryStore::create(state, reader);
    bk_library_tree_.set_model(bk_library_model_);
    bk_library_model_->appendColumnsTo(bk_library_tree_);
    bk_library_tree_.set_fixed_height_mode(true);
    bk_library_model_->set_sort_column(static_cast<LibraryStore::ColRec&>(*bk_library_model_->columns).id, Gtk::SortType::SORT_DESCENDING);
    bk_library_tree_.signal_row_activated().connect([this] (const Gtk::TreeModel::Path &path, Gtk::TreeViewColumn*) -> void {
        open_info_dialog_(bk_library_model_->member(path).id);
    });
    swins_.emplace_back();
    swins_.back().add(bk_library_tree_);
    swins_.back().set_policy(Gtk::POLICY_NEVER, Gtk::POLICY_AUTOMATIC);
    fields_["library_tab"].second.set_text("Library");
    tabs.append_page(swins_.back(), fields_["library_tab"].second);
    

    add(tabs);

    refresh(state);

    show_all();
}

InfoWindow::InfoWindow(std::shared_ptr<const State> state, std::shared_ptr<Gtk::Window> main_window, unsigned long book_id)
    : reader{0}, book{book_id}
{
    initWindow(*main_window);

    set_title("Book details (" + std::to_string(book) + ")");

    auto &grid_status = new_grid();

    data_append(grid_status, "id", "ID");
    data_append(grid_status, "position", "Position");
    data_append(grid_status, "market", "Market");
    data_append(grid_status, "price", "Price");
    data_append(grid_status, "quality", "Quality");
    data_append(grid_status, "created", "Period written");
    data_append(grid_status, "age", "Age");
    data_append(grid_status, "revenue", "Revenue (lifetime)");
    data_append(grid_status, "revenueLast", "Revenue (current)");
    data_append(grid_status, "sales", "Copies sold (lifetime)");
    data_append(grid_status, "salesLast", "Copies sold (current)");
    data_append(grid_status, "pirated", "Copies pirated (lifetime)");
    data_append(grid_status, "piratedLast", "Copies pirated (current)");
    data_append(grid_status, "copies", "Copies in world");
    data_append(grid_status, "author", "Author ID");
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
        updateValue("booksPurchased", r.library_purchased);
        updateValue("booksPirated", r.library_pirated);
        updateValue("booksNew", r.new_books.size());
        updateValue("booksNewPurchased", r.library_purchased_new);
        updateValue("booksNewPirated", r.library_pirated_new);
        updateValue("booksWritten", r.wrote.size());
        updateValue("bookLast", r.wrote.empty()
                ? "(never written)"
                : std::to_string(state->t - state->books.at(*r.wrote.crbegin()).created));
        updateValue("numFriends", r.friends.size());

        updateValue("library_tab", "Library (" + std::to_string(r.library.size() - r.wrote.size()) + ")");
        updateValue("authored_books_tab", "Authored Books (" + std::to_string(r.wrote.size()) + ")");

#define UPDATE_LIN(PREFIX, VAR) \
        updateValue(PREFIX + std::string("n"), VAR.n()); \
        updateValue(PREFIX + std::string("s2"), VAR.s2()); \
        updateMatrix(PREFIX + std::string("V"), VAR.s2() * VAR.Vinv().fullPivHouseholderQr().inverse(), true); \
        for (size_t i = 0; i < VAR.K(); i++) \
            updateValue(PREFIX + std::to_string(i), VAR.beta()[i]);
#define UPDATE_LIN_RB(PREFIX, BELIEF) UPDATE_LIN(PREFIX, r.BELIEF)

        UPDATE_LIN_RB("q_", quality);
        UPDATE_LIN_RB("p_", profit);
        UPDATE_LIN_RB("d_", demand);
        UPDATE_LIN_RB("pe_", profitExtrap());


        updateValue("_p_draws", r.profit.draw_rejection_success);
        updateValue("_p_discards", r.profit.draw_rejection_discards);
        updateValue("_p_drawtype", r.profit.last_draw_mode == BayesianLinearRestricted::DrawMode::Gibbs ? "gibbs" : "rejection");
        updateValue("_pe_draws", r.profitExtrap().draw_rejection_success);
        updateValue("_pe_discards", r.profitExtrap().draw_rejection_discards);
        updateValue("_pe_drawtype", r.profitExtrap().last_draw_mode == BayesianLinearRestricted::DrawMode::Gibbs ? "gibbs" : "rejection");
        updateValue("_d_draws", r.demand.draw_rejection_success);
        updateValue("_d_discards", r.demand.draw_rejection_discards);
        updateValue("_d_drawtype", r.demand.last_draw_mode == BayesianLinearRestricted::DrawMode::Gibbs ? "gibbs" : "rejection");

        for (unsigned long a : Reader::profit_stream_ages) {
            std::string code_prefix = "ps" + std::to_string(a) + "_";

            if (r.profit_stream.count(a)) {
                auto &psi = r.profit_stream.at(a);
                UPDATE_LIN(code_prefix, psi);
            }
            else {
                updateValue(code_prefix + "n", "0");
                updateValue(code_prefix + "s2", "");
                clearMatrix(code_prefix + "V");
                for (size_t j = 0; j < a; j++)
                    updateValue(code_prefix + std::to_string(j), "");
            }
        }

        // Update the books trees
        if (not initial_refresh_) {
            // Preserve the current sort column/order, if any
            int sort_col;
            Gtk::SortType sort_order;
            bool resort = bk_authored_model_->get_sort_column_id(sort_col, sort_order);
            bk_authored_model_ = BookStore::create(state, reader);
            if (resort) bk_authored_model_->set_sort_column(sort_col, sort_order);
            bk_authored_tree_.set_model(bk_authored_model_);

            resort = bk_library_model_->get_sort_column_id(sort_col, sort_order);
            bk_library_model_ = LibraryStore::create(state, reader);
            if (resort) bk_library_model_->set_sort_column(sort_col, sort_order);
            bk_library_tree_.set_model(bk_library_model_);
        }
    }
    else {
        if (state->books.count(book)) {
            set_title("Book details (" + std::to_string(book) + ")");
            auto &b = state->books.at(book);
            updateValue("id", b.id);
            updateValue("author", b.author);
            updateValue("position", GUI::pos_to_string(b.position));
            updateValue("market", b.market() ? "yes" : "no");
            updateValue("price", b.price);
            updateValue("quality", b.quality);
            updateValue("created", b.created);
            updateValue("age", state->t - b.created);
            updateValue("revenue", b.revenue_lifetime);
            updateValue("revenueLast", b.revenue);
            updateValue("sales", b.sales_lifetime);
            updateValue("salesLast", b.sales);
            updateValue("pirated", b.pirated_lifetime);
            updateValue("piratedLast", b.pirated);
            updateValue("copies", b.copies_lifetime());
        }
        else {
            // If the user navigates back in time to a period where the book doesn't exist, set
            // everything to N/A and change the title.
            set_title("Book details (" + std::to_string(book) + "): not yet created!");
            for (auto &key : {"id", "author", "position", "market", "price", "quality", "age", "created",
                    "revenue", "revenueLast", "sales", "salesLast", "pirated", "piratedLast", "copies"}) {
                updateValue(key, "N/A");
            }
        }
    }

    if (initial_refresh_) initial_refresh_ = false;
}

void InfoWindow::updateValue(const std::string &code, const std::string &val) {
    fields_[code].second.set_text(val);
}

void InfoWindow::updateMatrix(const std::string &code, const Ref<const MatrixXd> &m, bool lower_triangle) {
    size_t pos = 0;
    auto &labels = matrix_[code];
    for (int i = 0; i < m.rows(); i++) for (int j = 0; j < m.cols(); j++) {
        if (pos >= labels.size()) return;
        std::ostringstream sstr;
        if (not lower_triangle or j <= i) sstr << m(i,j);
        labels[pos]->set_markup(sstr.str());
        pos++;
    }
}

void InfoWindow::clearMatrix(const std::string &code, const std::string &value) {
    for (auto &label : matrix_[code]) label->set_markup(value);
}

} }
