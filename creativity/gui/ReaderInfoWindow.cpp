#include "creativity/gui/ReaderInfoWindow.hpp"
#include "creativity/gui/GUI.hpp"
#include "creativity/gui/BookStore.hpp"
#include "creativity/gui/LibraryStore.hpp"
#include "creativity/Reader.hpp"
#include "creativity/belief/Demand.hpp"
#include "creativity/belief/Profit.hpp"
#include "creativity/belief/ProfitStream.hpp"
#include "creativity/state/State.hpp"
#include <Eigen/Core>
#include <Eigen/QR>
#include <glibmm/signalproxy.h>
#include <gtkmm/enums.h>
#include <gtkmm/label.h>
#include <gtkmm/notebook.h>
#include <gtkmm/scrolledwindow.h>
#include <gtkmm/treemodel.h>
#include <gtkmm/treemodelcolumn.h>
#include <gtkmm/treeview.h>
#include <gtkmm/treeviewcolumn.h> // IWYU pragma: keep
#include <cstddef>
#include <iterator>
#include <list>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>


namespace creativity { namespace gui {

using namespace eris;
using namespace eris::belief;
using namespace creativity::state;

ReaderInfoWindow::ReaderInfoWindow(
        std::shared_ptr<const state::State> state,
        std::shared_ptr<Gtk::Window> main_window,
        eris_id_t reader,
        std::function<void(eris::eris_id_t)> open_info_dialog
        )
    : InfoWindow(std::move(main_window), reader, "Reader/author details (" + std::to_string(reader) + ")"),
    open_info_dialog_{std::move(open_info_dialog)}
{
    if (state->readers.count(reader) == 0)
        throw std::out_of_range("ReaderInfoWindow() called with invalid reader id");

    auto rdr = state->readers.at(reader);

    nbs_.emplace_back();
    Gtk::Notebook &tabs = nbs_.back();

    auto &grid_status = new_tab_grid(tabs, "Status");
    data_append(grid_status, "id", "ID");
    data_append(grid_status, "position", "Position");
    data_append(grid_status, "utility", "Utility");
    data_append(grid_status, "uLife", "Lifetime utility");
    data_append(grid_status, "creationScale", "Creation scale");
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

    const std::string beta(u8"<i>β</i>");

    auto &grid_profit = new_tab_grid(beliefs, "Profit");
    labels_append(grid_profit, "Dependent variable", "<i>lifetimeProfit</i>");
    data_append(grid_profit, "p_n", "n");
    data_append(grid_profit, "p_s2", "s<sup>2</sup>");
    std::vector<std::string> p_vars(rdr.profit.names());
    for (size_t i = 0; i < p_vars.size(); i++)
        data_append(grid_profit, "p_" + std::to_string(i), beta + "[" + p_vars[i] + "]");
    data_append(grid_profit, "_p_draws", "# successful draws");
    data_append(grid_profit, "_p_discards", "# discarded draws");
    data_append(grid_profit, "_p_drawtype", "Draw method");
    matrix_at(grid_profit, "p_s2V", "s<sup>2</sup><b>V</b>", 2, 2, p_vars.size(), p_vars.size());

    auto &grid_demand = new_tab_grid(beliefs, "Demand");
    labels_append(grid_demand, "Dependent variable", "<i>quantityDemanded</i>");
    data_append(grid_demand, "d_n", "n");
    data_append(grid_demand, "d_s2", "s<sup>2</sup>");
    std::vector<std::string> d_vars(rdr.demand.names());
    for (size_t i = 0; i < d_vars.size(); i++)
        data_append(grid_demand, "d_" + std::to_string(i), beta + "[" + d_vars[i] + "]");
    data_append(grid_demand, "_d_draws", "# successful draws");
    data_append(grid_demand, "_d_discards", "# discarded draws");
    data_append(grid_demand, "_d_drawtype", "Draw method");
    matrix_at(grid_demand, "d_s2V", "s<sup>2</sup><b>V</b>", 2, 2, d_vars.size(), d_vars.size());
    comment_append(grid_demand, "<i>NB: This regression is for single-period demand.</i>", p_vars.size() + 3);

    nbs_.emplace_back();
    Gtk::Notebook &pstream_tabs = nbs_.back();
    beliefs.append_page(pstream_tabs, "Profit Stream");


    for (unsigned long a : Reader::profit_stream_ages) {
        auto &grid_pstream = new_tab_grid(pstream_tabs, u8"age ≥ " + std::to_string(a));
        labels_append(grid_pstream, "Dependent variable", "<i>profitRemaining</i>");
        std::string code_prefix = "ps" + std::to_string(a) + "_";

        data_append(grid_pstream, code_prefix + "n", "n");
        data_append(grid_pstream, code_prefix + "s2", "s<sup>2</sup>");
        for (unsigned long j = 0; j < a; j++) {
            data_append(grid_pstream, code_prefix + std::to_string(j), beta + u8"[π<sub>" + std::to_string(j) + "</sub>]");
        }
        matrix_at(grid_pstream, code_prefix + "s2V", "s<sup>2</sup><b>V</b>", 2, 2, a, a);
    }

    auto &grid_pextrap = new_tab_grid(beliefs, "Profit (extrap.)");
    labels_append(grid_pextrap, "Dependent variable", "<i>lifetimeProfit</i>");
    data_append(grid_pextrap, "pe_n", "n");
    data_append(grid_pextrap, "pe_s2", "s<sup>2</sup>");
    for (size_t i = 0; i < p_vars.size(); i++)
        data_append(grid_pextrap, "pe_" + std::to_string(i), beta + "[" + p_vars[i] + "]");
    data_append(grid_pextrap, "_pe_draws", "# successful draws");
    data_append(grid_pextrap, "_pe_discards", "# discarded draws");
    data_append(grid_pextrap, "_pe_drawtype", "Draw method");
    matrix_at(grid_pextrap, "pe_s2V", "s<sup>2</sup><b>V</b>", 2, 2, p_vars.size(), p_vars.size());
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

void ReaderInfoWindow::refresh(std::shared_ptr<const State> state) {
    if (state->t == t_) return;

    auto &r = state->readers.at(id);

    updateValue("id", r.id);
    updateValue("position", GUI::pos_to_string(r.position));
    updateValue("utility", r.u);
    updateValue("uLife", r.u_lifetime);
    updateValue("creationScale", r.creation_scale);
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
    updateMatrix(PREFIX + std::string("s2V"), VAR.s2() * VAR.Vinv().fullPivHouseholderQr().inverse(), true); \
    for (size_t i = 0; i < VAR.K(); i++) \
        updateValue(PREFIX + std::to_string(i), VAR.beta()[i]);
#define UPDATE_LIN_RB(PREFIX, BELIEF) UPDATE_LIN(PREFIX, r.BELIEF)

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
            clearMatrix(code_prefix + "s2V");
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
        bk_authored_model_ = BookStore::create(state, id);
        if (resort) bk_authored_model_->set_sort_column(sort_col, sort_order);
        bk_authored_tree_.set_model(bk_authored_model_);

        resort = bk_library_model_->get_sort_column_id(sort_col, sort_order);
        bk_library_model_ = LibraryStore::create(state, id);
        if (resort) bk_library_model_->set_sort_column(sort_col, sort_order);
        bk_library_tree_.set_model(bk_library_model_);
    }

    if (initial_refresh_) initial_refresh_ = false;
}

}}
