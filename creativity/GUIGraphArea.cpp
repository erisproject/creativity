#include "creativity/GUIGraphArea.hpp"
#include "creativity/GUI.hpp"
#include "creativity/Reader.hpp"
#include "creativity/Book.hpp"
#include <eris/Random.hpp>
#include <cmath>

using namespace eris;

namespace creativity {

GUIGraphArea::GUIGraphArea(const double &top, const double &right, const double &bottom, const double &left,
        std::shared_ptr<Simulation> sim, GUI &gui) :
    sim_{sim}, gui_{gui}, bounds_{top,right,bottom,left}
{}

Cairo::Matrix GUIGraphArea::graph_to_canvas() const {
    Gtk::Allocation allocation = get_allocation();
    const int width = allocation.get_width();
    const int height = allocation.get_height();

    // The second of these is typically negative, which is correct because "up" in the graph
    // translates to a lower coordinate on the screen (since positive y coordinates are down the
    // screen).
    const double gwidth = bounds_.right - bounds_.left,
          gheight = bounds_.bottom - bounds_.top;

    // Build a transformation that converts from positions to canvas coordinates
    auto trans = Cairo::identity_matrix();
    trans.scale(width / gwidth, height / gheight);
    trans.translate(-bounds_.left, -bounds_.top);

    return trans;
}

bool GUIGraphArea::on_draw(const Cairo::RefPtr<Cairo::Context> &cr) {
    ERIS_DBG("Starting draw");
    Gtk::Allocation allocation = get_allocation();
    const int width = allocation.get_width();
    const int height = allocation.get_height();

    // Lock out the simulation until we finish redrawing
    auto lock = sim_->runLock();
    ERIS_DBG("Got sim lock");

    cr->save();

    // Paint the background white
    cr->set_source_rgb(1,1,1);
    cr->paint();

    auto trans = graph_to_canvas();

    // Add axes
    double zero_x = 0, zero_y = 0;
    trans.transform_point(zero_x, zero_y);
    cr->set_line_width(1.0);
    cr->set_source_rgb(0,0,0);
    cr->move_to(zero_x, 0);
    cr->line_to(zero_x, height);
    cr->stroke();
    cr->move_to(0, zero_y);
    cr->line_to(width, zero_y);
    cr->stroke();

    // Tick marks
    cr->set_line_width(1.0);
    const double tick_max_x = std::max(std::abs(bounds_.right), std::abs(bounds_.left));
    int tick_num = 0;
    for (double gx = tick_space; gx <= tick_max_x; gx += tick_space) {
        const double curr_tick_size = (++tick_num % tick_big) ? tick_size : 3*tick_size;
        double x = gx, y = 0;
        trans.transform_point(x, y);
        cr->move_to(x, y - curr_tick_size/2.0);
        cr->rel_line_to(0, curr_tick_size);
        x = -gx; y = 0;
        trans.transform_point(x, y);
        cr->move_to(x, y - curr_tick_size/2.0);
        cr->rel_line_to(0, curr_tick_size);
    }
    const double tick_max_y = std::max(std::abs(bounds_.top), std::abs(bounds_.bottom));
    tick_num = 0;
    for (double gy = tick_space; gy <= tick_max_y; gy += tick_space) {
        const double curr_tick_size = (++tick_num % tick_big) ? tick_size : 3*tick_size;
        double x = 0, y = gy;
        trans.transform_point(x, y);
        cr->move_to(x - curr_tick_size/2.0, y);
        cr->rel_line_to(curr_tick_size, 0);
        x = 0; y = -gy;
        trans.transform_point(x, y);
        cr->move_to(x - curr_tick_size/2.0, y);
        cr->rel_line_to(curr_tick_size, 0);
    }
    cr->stroke();

    SharedMember<Reader> token_reader; // Need to keep a reader for wrapping

    ERIS_DBG("Drawing readers");
    cr->set_source_rgba(1, 0.55, 0, 0.5); // Lines from readers to purchased books
    for (auto &r : sim_->agents<Reader>()) {
        if (!token_reader) { token_reader = r; }
        double rx = r->position()[0], ry = r->position()[1];
        drawPoint(cr, trans, rx, ry, PointType::X, 1, 0, 0, 1, r);

        const double radius1 = 0.05 * (r->u() - 1000);
        if (radius1 > 0) { // Utility > 1000 means some utility gain from books
            drawCircle(cr, trans, rx, ry, radius1, .133, .545, .133, 0.5);
        }
        const double radius2 = 0.1 * r->newBooks().size();
        if (radius2 > 0) { // Bought some books
            drawCircle(cr, trans, rx, ry, radius2, 1, .55, 0, 0.5);
            for (auto &book_id : r->newBooks()) {
                auto b = sim_->good<Book>(book_id);
                drawWrappingLine(cr, trans, r, b);
            }
            cr->stroke();
        }
        // NB: rx, ry may be translated now
    }

    cr->set_source_rgba(0.5, 0.2, 0.5, 0.5); // Colour for lines from books to their author
    for (auto &b : sim_->goods<Book>()) {
        // Give new books a larger cross: a brand new book gets a point 3 times as large; this
        // scaling decreases linearly to the regular size at a 10-period old book.
        const double scale = std::max(1.0, 3.0 - 0.3*b->age());

        // Draw book:
        double br = 0, bg = .4, bb = 1; // Blue with a hint of green
        if (not b->hasMarket())
            br = bg = bb = 0;
        drawPoint(cr, trans, b->position()[0], b->position()[1], PointType::CROSS, br, bg, bb, 1, token_reader, scale);

        drawWrappingLine(cr, trans, b->author(), b);

        cr->stroke();
    }

    cr->restore();

    auto it = gui_.info_windows_.begin();
    while (it != gui_.info_windows_.end()) {
        if (not it->second.get_visible()) {
            it = gui_.info_windows_.erase(it);
        }
        else {
            it->second.refresh();
            it++;
        }
    }

    gui_.queueEvent(GUI::Event::Type::redraw);
    return true;
}

void GUIGraphArea::drawWrappingLine(const Cairo::RefPtr<Cairo::Context> &cr, const Cairo::Matrix &trans, const Reader &r, const Book &b) {
    cr->save();
    cr->transform(trans);
    const auto &rp = r.position();
    const double x_span = r.wrapUpperBound()[0] - r.wrapLowerBound()[0];
    const double y_span = r.wrapUpperBound()[1] - r.wrapLowerBound()[1];
    auto v = r.vectorTo(b);
    // There are nine virtual points the author can take; draw a line from each of them (at most
    // three of these lines will actually show up)
    for (double r_x : {rp[0], rp[0] - x_span, rp[0] + x_span}) {
        for (double r_y : {rp[1], rp[1] - y_span, rp[1] + y_span}) {
            cr->move_to(r_x, r_y);
            cr->rel_line_to(v[0], v[1]);
        }
    }
    cr->restore(); // Undo transformation
}

#define RGBA red, green, blue, alpha
void GUIGraphArea::drawPoint(
        const Cairo::RefPtr<Cairo::Context> &cr, const Cairo::Matrix &trans,
        double x, double y,
        const PointType &type, double red, double green, double blue, double alpha,
        const SharedMember<Reader> &r, double scale) {

    // The radius of the point
    const double pt_radius = point_size * scale;

    // Before we transform x and y, consider whether there are wrapped versions of the point we need
    // to worry about; if so, draw them first.
    if (r) {
        SharedMember<Reader> null_reader;
        double dim_x = r->wrapUpperBound()[0] - r->wrapLowerBound()[0],
               dim_y = r->wrapUpperBound()[1] - r->wrapLowerBound()[1];
        bool wrap_right = x + pt_radius > r->wrapUpperBound()[0],
             wrap_left = x - pt_radius < r->wrapLowerBound()[0],
             wrap_up = y + pt_radius > r->wrapUpperBound()[1],
             wrap_down = y - pt_radius < r->wrapLowerBound()[1];

        // There are 8 wrapping possibilities to consider; the first 4 are corners, the last 4 are
        // single-dimension edge wraps.  Corners are just like double-edged wraps, but *also* have
        // to wrap to the opposite corner.
        if (wrap_left and wrap_up) // Upper left
            drawPoint(cr, trans, x+dim_x, y-dim_y, type, RGBA, null_reader, scale); // Lower right
        else if (wrap_left and wrap_down) // Lower left
            drawPoint(cr, trans, x+dim_x, y+dim_y, type, RGBA, null_reader, scale); // Upper right
        else if (wrap_right and wrap_up) // Upper right
            drawPoint(cr, trans, x-dim_x, y-dim_y, type, RGBA, null_reader, scale); // Lower left
        else if (wrap_right and wrap_down) // Lower right
            drawPoint(cr, trans, x-dim_x, y+dim_y, type, RGBA, null_reader, scale); // Upper left

        // Now the edge mirroring
        if (wrap_left)
            drawPoint(cr, trans, x+dim_x, y, type, RGBA, null_reader, scale);
        else if (wrap_right)
            drawPoint(cr, trans, x-dim_x, y, type, RGBA, null_reader, scale);
        
        if (wrap_up)
            drawPoint(cr, trans, x, y-dim_y, type, RGBA, null_reader, scale);
        else if (wrap_down)
            drawPoint(cr, trans, x, y+dim_y, type, RGBA, null_reader, scale);
    }

    // Get the user-space coordinates of the point
    trans.transform_point(x, y);

    cr->save();
    cr->set_line_width(2.0);

    switch (type) {
        case PointType::X:
            {
                // We're drawing 45-degree lines, so get the linear distance on each dimension:
                double edge = pt_radius * sqrt(2.0);

                cr->set_source_rgba(RGBA);

                cr->move_to(x-0.5*edge, y-0.5*edge); // Upper left
                cr->rel_line_to(edge, edge); // ... line to bottom right

                cr->rel_move_to(-edge, 0); // Lower left
                cr->rel_line_to(edge, -edge); // ... to upper right
                break;
            }
        case PointType::SQUARE:
            {
                // The SQUARE vertices are the the same as the X type, above
                double edge = pt_radius * sqrt(2.0);

                cr->set_source_rgba(RGBA);
                cr->rectangle(x-0.5*edge, y-0.5*edge, edge, edge);
                break;
            }
        case PointType::CROSS:
            {
                cr->set_source_rgba(RGBA);
                // Horizontal line:
                cr->move_to(x - pt_radius, y);
                cr->rel_line_to(2*pt_radius, 0);

                // Verical line:
                cr->move_to(x, y - pt_radius);
                cr->rel_line_to(0, 2*pt_radius);
                break;
            }
    }

    cr->stroke();
    cr->restore();
}

void GUIGraphArea::drawCircle(const Cairo::RefPtr<Cairo::Context> &cr, const Cairo::Matrix &trans,
        double cx, double cy, double r, double red, double green, double blue, double alpha) {
    // Apply the full transformation: then we can just draw the circle in *graph* space which,
    // post-transformation, will be exactly the desired oval.
    cr->save();
    cr->transform(trans);

    // Set up the circle in the transformed (graph) coordinates, but *don't* draw it yet (if we do,
    // the line width will also be transformed, which we don't want--it would be too thin on the
    // more compressed axis).
    cr->arc(cx, cy, r, 0.0, 2*M_PI);

    // Draw a radial line at a quasi-random angle (based on the position and radius values) to the
    // center of the circle
    //std::uniform_real_distribution<double> rand_angle(0, 2*M_PI);
    //double angle = rand_angle(Random::rng());
    double angle = std::fmod(1000000 * cx * cy * r, 2*M_PI);
    cr->move_to(cx, cy);
    cr->rel_line_to(r*std::sin(angle), r*std::cos(angle));

    // Undo the transformation, then change the color and line width and actually draw the circle
    cr->restore();
    cr->save();
    cr->set_line_width(2.0);
    cr->set_source_rgba(RGBA);

    cr->stroke();
    cr->restore();
}

}

