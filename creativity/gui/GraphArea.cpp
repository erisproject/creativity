#include "creativity/gui/GraphArea.hpp"
#include "creativity/gui/GUI.hpp"
#include "creativity/Creativity.hpp"
#include "creativity/Reader.hpp"
#include "creativity/Book.hpp"
#include <cmath>

using namespace eris;
using namespace creativity::state;

namespace creativity { namespace gui {

GraphArea::GraphArea(GUI &gui) :
    gui_{gui}, wpb_({0.0,0.0})
{}

Cairo::Matrix GraphArea::graph_to_canvas() const {
    Gtk::Allocation allocation = get_allocation();
    const int width = allocation.get_width();
    const int height = allocation.get_height();

    const double boundary = gui_.creativity_->boundary();

    // The second of these is negative, which is correct because "up" in the graph translates to a
    // lower coordinate on the screen (since positive y coordinates are down the screen).
    const double gwidth = 2*boundary, gheight = -2*boundary;

    // Build a transformation that converts from positions to canvas coordinates
    auto trans = Cairo::identity_matrix();
    trans.scale(width / gwidth, height / gheight);
    trans.translate(boundary, -boundary);

    return trans;
}

bool GraphArea::on_draw(const Cairo::RefPtr<Cairo::Context> &cr_grapharea) {
    Gtk::Allocation allocation = get_allocation();
    const int width = allocation.get_width();
    const int height = allocation.get_height();
    const double boundary = gui_.creativity_->boundary();

    if (gui_.state_num_ == 0) {
        // No states at all: just draw a blank screen
        cr_grapharea->set_source(colours.background);
        cr_grapharea->paint();
        return true;
    }

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

        std::shared_ptr<const State> state;

        {
            auto st = gui_.creativity_->storage();
            if (st.first->empty())
                state = std::make_shared<State>();
            else
                state = (*st.first)[gui_.state_curr_];
        }

        cr->save();

        // Paint the background white
        cr->set_source(colours.background);
        cr->paint();

        auto trans = graph_to_canvas();

        // Things get draw here in order from least to most important (since later things are drawn
        // on top of earlier things).  This means we have to loop over the same vector a few times,
        // unfortunately.

        // Draw lines from readers to newly purchased books.
        if (not invisible(colours.reading)) {
            cr->set_source(colours.reading);
            for (auto &rpair : state->readers) {
                auto &r = rpair.second;

                for (auto &book_id : r.new_books) {
                    auto &b = state->books.at(book_id);
                    drawWrappingLine(cr, trans, r.position, b.position);
                }
            }
            cr->stroke();
        }

        // Draw books.  On-market and off-market books have different colours, and newer books are
        // larger than older books.
        if (not invisible(colours.book_live) or not invisible(colours.book_dead)) {
            for (auto &bpair : state->books) {
                auto &b = bpair.second;
                // Give new books a larger cross: a brand new book gets a point 3 times as large; this
                // scaling decreases linearly to the regular size at a 10-period old book.
                const double scale = std::max(1.0, 3.0 - 0.3*b.age);

                // Draw book:
                drawPoint(cr, trans, b.position[0], b.position[1], PointType::CROSS,
                        b.market ? colours.book_live : colours.book_dead, scale);
            }
        }

        // Readers have utility circles, indicating how much about the base 1000 utility the reader
        // was in the period.
        if (not invisible(colours.utility)) {
            for (auto &rpair : state->readers) {
                auto &r = rpair.second;
                double rx = r.position[0], ry = r.position[1];
                const double radius1 = 0.05 * (r.u - 1000);
                if (radius1 > 0) { // Utility > 1000 means some utility gain from books
                    drawCircle(cr, trans, rx, ry, radius1, colours.utility);
                }
            }
        }

        // Lines from each book to its author
        if (not invisible(colours.author_live) or not invisible(colours.author_dead)) {
            for (auto &bpair : state->books) {
                auto &b = bpair.second;
                cr->set_source(b.market ? colours.author_live : colours.author_dead);
                drawWrappingLine(cr, trans, state->readers.at(b.author).position, b.position);
                cr->stroke();
            }
        }

        // Draw readers
        if (not invisible(colours.reader)) {
            for (auto &rpair : state->readers) {
                auto &r = rpair.second;
                double rx = r.position[0], ry = r.position[1];
                // Draw the reader
                drawPoint(cr, trans, rx, ry, PointType::X, colours.reader);

                // NB: rx, ry may be translated now
            }
        }

        // Add axes (last, so that they are on top)
        double zero_x = 0, zero_y = 0;
        trans.transform_point(zero_x, zero_y);

        // Round the transformed location to the nearest value; it'll be ever so slighly inaccurate,
        // but will give a double pixel line instead of a slightly blurry line (from partial pixel coverage).
        zero_x = std::round(zero_x);
        zero_y = std::round(zero_y);

        cr->set_line_width(2.0);
        cr->set_source(colours.axes);
        cr->move_to(zero_x, 0);
        cr->line_to(zero_x, height);
        cr->stroke();
        cr->move_to(0, zero_y);
        cr->line_to(width, zero_y);
        cr->stroke();

        // Tick marks
        cr->set_line_width(1.0);
        int tick_num = 0;
        for (double gx = tick_space; gx <= boundary; gx += tick_space) {
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
        tick_num = 0;
        for (double gy = tick_space; gy <= boundary; gy += tick_space) {
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

        cr->restore();
    }

    cr_grapharea->set_source(drawing_cache_[gui_.state_curr_], 0, 0);
    cr_grapharea->paint();

    return true;
}

void GraphArea::drawWrappingLine(const Cairo::RefPtr<Cairo::Context> &cr, const Cairo::Matrix &trans, const Position &from, const Position &to) {
    cr->save();
    cr->transform(trans);
    const double boundary = gui_.creativity_->boundary();
    const double x_span = 2*boundary;
    const double y_span = 2*boundary;
    if (not wpb_.wrapped(0)) {
        wpb_ = WrappedPositionalBase({0.0,0.0}, {-boundary, -boundary}, {boundary, boundary});
    }
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
        const PointType &type, const Colour &colour, double scale, bool virt) {

    // The radius of the point
    const double pt_radius = point_size * scale;
    const double boundary = gui_.creativity_->boundary();

    // Before we transform x and y, consider whether there are wrapped versions of the point we need
    // to worry about; if so, draw them first.
    if (not virt) {
        double dim_x = 2*boundary,
               dim_y = 2*boundary;
        bool wrap_right = x + pt_radius > boundary,
             wrap_left  = x - pt_radius < -boundary,
             wrap_up    = y + pt_radius > boundary,
             wrap_down  = y - pt_radius < -boundary;

        // There are 8 wrapping possibilities to consider; the first 4 are corners, the last 4 are
        // single-dimension edge wraps.  Corners are just like double-edged wraps, but *also* have
        // to wrap to the opposite corner.
        if (wrap_left and wrap_up) // Upper left
            drawPoint(cr, trans, x+dim_x, y-dim_y, type, colour, scale, true); // Lower right
        else if (wrap_left and wrap_down) // Lower left
            drawPoint(cr, trans, x+dim_x, y+dim_y, type, colour, scale, true); // Upper right
        else if (wrap_right and wrap_up) // Upper right
            drawPoint(cr, trans, x-dim_x, y-dim_y, type, colour, scale, true); // Lower left
        else if (wrap_right and wrap_down) // Lower right
            drawPoint(cr, trans, x-dim_x, y+dim_y, type, colour, scale, true); // Upper left

        // Now the edge mirroring
        if (wrap_left)
            drawPoint(cr, trans, x+dim_x, y, type, colour, scale, true);
        else if (wrap_right)
            drawPoint(cr, trans, x-dim_x, y, type, colour, scale, true);
        
        if (wrap_up)
            drawPoint(cr, trans, x, y-dim_y, type, colour, scale, true);
        else if (wrap_down)
            drawPoint(cr, trans, x, y+dim_y, type, colour, scale, true);
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

                cr->set_source(colour);

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

                cr->set_source(colour);
                cr->rectangle(x-0.5*edge, y-0.5*edge, edge, edge);
                break;
            }
        case PointType::CROSS:
            {
                cr->set_source(colour);
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
        double cx, double cy, double r, const Colour &colour) {
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
    cr->set_source(colour);

    cr->stroke();
    cr->restore();
}

} }

