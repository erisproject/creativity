#include "creativity/gui/GraphArea.hpp"
#include "creativity/gui/GUI.hpp"
#include "creativity/Reader.hpp"
#include "creativity/Book.hpp"
#include <cmath>

using namespace eris;
using namespace creativity::state;

namespace creativity { namespace gui {

GraphArea::GraphArea(const double &top, const double &right, const double &bottom, const double &left, GUI &gui) :
    gui_{gui}, wpb_({0.0,0.0}, {-BOUNDARY, -BOUNDARY}, {BOUNDARY, BOUNDARY}), bounds_{top,right,bottom,left}
{}

Cairo::Matrix GraphArea::graph_to_canvas() const {
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

bool GraphArea::on_draw(const Cairo::RefPtr<Cairo::Context> &cr_grapharea) {
    Gtk::Allocation allocation = get_allocation();
    const int width = allocation.get_width();
    const int height = allocation.get_height();

    if (drawing_cache_width_ != width || drawing_cache_height_ != height) {
        drawing_cache_.clear();
        drawing_cache_width_ = width;
        drawing_cache_height_ = height;
    }

    if (drawing_cache_.size() < gui_.state_curr_ + 1) {
        drawing_cache_.resize(gui_.state_num_);
    }
    if (not drawing_cache_[gui_.state_curr_]) {
        auto cr = Cairo::Context::create(
                drawing_cache_[gui_.state_curr_] = Cairo::ImageSurface::create(Cairo::FORMAT_RGB24, width, height)
                );

        std::shared_ptr<State> state;

        {
            // Lock the state vector while we grab the state
            auto lock = gui_.stateLock();
            if (gui_.states_.empty())
                state = std::make_shared<State>();
            else
                state = gui_.states_[gui_.state_curr_];
        }

        cr->save();

        // Paint the background white
        cr->set_source_rgb(1,1,1);
        cr->paint();

        auto trans = graph_to_canvas();

        // Add axes
        double zero_x = 0, zero_y = 0;
        trans.transform_point(zero_x, zero_y);
        // Round the transformed location to the nearest .5 value; it'll be ever so slighly inaccurate,
        // but will give a single pixel line instead of a double-width, blurry line that results from
        // positioning at an integer.
        zero_x = std::round(zero_x + 0.5) - 0.5;
        zero_y = std::round(zero_y + 0.5) - 0.5;

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
            x = std::round(x + 0.5) - 0.5; // Round to nearest 0.5 for crisper lines
            cr->move_to(x, y - curr_tick_size/2.0);
            cr->rel_line_to(0, curr_tick_size);
            x = -gx; y = 0;
            trans.transform_point(x, y);
            x = std::round(x + 0.5) - 0.5; // Round for crisper lines
            cr->move_to(x, y - curr_tick_size/2.0);
            cr->rel_line_to(0, curr_tick_size);
        }
        const double tick_max_y = std::max(std::abs(bounds_.top), std::abs(bounds_.bottom));
        tick_num = 0;
        for (double gy = tick_space; gy <= tick_max_y; gy += tick_space) {
            const double curr_tick_size = (++tick_num % tick_big) ? tick_size : 3*tick_size;
            double x = 0, y = gy;
            trans.transform_point(x, y);
            y = std::round(y + 0.5) - 0.5; // Round for crisper lines
            cr->move_to(x - curr_tick_size/2.0, y);
            cr->rel_line_to(curr_tick_size, 0);
            x = 0; y = -gy;
            trans.transform_point(x, y);
            y = std::round(y + 0.5) - 0.5; // Round for crisper lines
            cr->move_to(x - curr_tick_size/2.0, y);
            cr->rel_line_to(curr_tick_size, 0);
        }
        cr->stroke();

        //SharedMember<Reader> token_reader; // Need to keep a reader for wrapping

        cr->set_source_rgba(1, 0.55, 0, 0.5); // Lines from readers to purchased books
        for (auto &rpair : state->readers) {
            auto &r = rpair.second;
            //if (!token_reader) { token_reader = r; }
            double rx = r.position[0], ry = r.position[1];
            drawPoint(cr, trans, rx, ry, PointType::X, 1, 0, 0, 1);

            const double radius1 = 0.05 * (r.u - 1000);
            if (radius1 > 0) { // Utility > 1000 means some utility gain from books
                drawCircle(cr, trans, rx, ry, radius1, .133, .545, .133, 0.5);
            }
            const double radius2 = 0.1 * r.newBooks.size();
            if (radius2 > 0) { // Bought some books
                drawCircle(cr, trans, rx, ry, radius2, 1, .55, 0, 0.5);
                for (auto &book_id : r.newBooks) {
                    auto &b = state->books.at(book_id);
                    drawWrappingLine(cr, trans, r.position, b.position);
                }
                cr->stroke();
            }
            // NB: rx, ry may be translated now
        }

        cr->set_source_rgba(0.5, 0.2, 0.5, 0.5); // Colour for lines from books to their author
        for (auto &bpair : state->books) {
            auto &b = bpair.second;
            // Give new books a larger cross: a brand new book gets a point 3 times as large; this
            // scaling decreases linearly to the regular size at a 10-period old book.
            const double scale = std::max(1.0, 3.0 - 0.3*b.age);

            // Draw book:
            double br = 0, bg = .4, bb = 1; // Blue with a hint of green
            if (not b.market)
                br = bg = bb = 0;
            drawPoint(cr, trans, b.position[0], b.position[1], PointType::CROSS, br, bg, bb, 1, scale);

            drawWrappingLine(cr, trans, state->readers.at(b.author).position, b.position);

            cr->stroke();
        }

        cr->restore();
    }

    cr_grapharea->set_source(drawing_cache_[gui_.state_curr_], 0, 0);
    cr_grapharea->paint();

    return true;
}

void GraphArea::drawWrappingLine(const Cairo::RefPtr<Cairo::Context> &cr, const Cairo::Matrix &trans, const Position &from, const Position &to) {
    cr->save();
    cr->transform(trans);
    const double x_span = 2*BOUNDARY;
    const double y_span = 2*BOUNDARY;
    wpb_.moveTo(from);
    auto v = wpb_.vectorTo(to);
    // There are nine virtual points the author can take; draw a line from each of them (at most
    // three of these lines will actually show up)
    for (double r_x : {from[0], from[0] - x_span, from[0] + x_span}) {
        for (double r_y : {from[1], from[1] - y_span, from[1] + y_span}) {
            cr->move_to(r_x, r_y);
            cr->rel_line_to(v[0], v[1]);
        }
    }
    cr->restore(); // Undo transformation
}

#define RGBA red, green, blue, alpha
void GraphArea::drawPoint(
        const Cairo::RefPtr<Cairo::Context> &cr, const Cairo::Matrix &trans, double x, double y,
        const PointType &type, double red, double green, double blue, double alpha, double scale, bool virt) {

    // The radius of the point
    const double pt_radius = point_size * scale;

    // Before we transform x and y, consider whether there are wrapped versions of the point we need
    // to worry about; if so, draw them first.
    if (not virt) {
        double dim_x = 2*BOUNDARY,
               dim_y = 2*BOUNDARY;
        bool wrap_right = x + pt_radius > BOUNDARY,
             wrap_left  = x - pt_radius < -BOUNDARY,
             wrap_up    = y + pt_radius > BOUNDARY,
             wrap_down  = y - pt_radius < -BOUNDARY;

        // There are 8 wrapping possibilities to consider; the first 4 are corners, the last 4 are
        // single-dimension edge wraps.  Corners are just like double-edged wraps, but *also* have
        // to wrap to the opposite corner.
        if (wrap_left and wrap_up) // Upper left
            drawPoint(cr, trans, x+dim_x, y-dim_y, type, RGBA, scale, true); // Lower right
        else if (wrap_left and wrap_down) // Lower left
            drawPoint(cr, trans, x+dim_x, y+dim_y, type, RGBA, scale, true); // Upper right
        else if (wrap_right and wrap_up) // Upper right
            drawPoint(cr, trans, x-dim_x, y-dim_y, type, RGBA, scale, true); // Lower left
        else if (wrap_right and wrap_down) // Lower right
            drawPoint(cr, trans, x-dim_x, y+dim_y, type, RGBA, scale, true); // Upper left

        // Now the edge mirroring
        if (wrap_left)
            drawPoint(cr, trans, x+dim_x, y, type, RGBA, scale, true);
        else if (wrap_right)
            drawPoint(cr, trans, x-dim_x, y, type, RGBA, scale, true);
        
        if (wrap_up)
            drawPoint(cr, trans, x, y-dim_y, type, RGBA, scale, true);
        else if (wrap_down)
            drawPoint(cr, trans, x, y+dim_y, type, RGBA, scale, true);
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

void GraphArea::drawCircle(const Cairo::RefPtr<Cairo::Context> &cr, const Cairo::Matrix &trans,
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
    // center of the circle.  The quasi-random calculation means the same circle gets the same angle
    // when redrawn.
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

} }

