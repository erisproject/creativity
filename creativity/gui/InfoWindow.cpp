#include "creativity/gui/InfoWindow.hpp"
#include "creativity/gui/GUI.hpp"
#include <iomanip>
#include <gtkmm/box.h>
#include <gtkmm/treemodelcolumn.h>
#include <Eigen/Core>

using namespace eris;
using namespace creativity::state;
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


InfoWindow::InfoWindow(std::shared_ptr<Gtk::Window> main_window, eris_id_t id, const std::string &title)
    : id{id}
{
    set_transient_for(*main_window);
//    set_type_hint(Gdk::WINDOW_TYPE_HINT_UTILITY);

    set_resizable(true);
    property_destroy_with_parent() = true;
    set_deletable(true);

    set_title(title);
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
