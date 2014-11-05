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

    const double half_width = 0.5 * allocation.get_width(), half_height = 0.5 * allocation.get_height();
    const double boundary = gui_.creativity_->boundary();
    // The second of these is negative because higher y canvas values are down the screen, which are
    // *lower* graph values.
    return Cairo::Matrix(half_width / boundary, 0, 0, -half_height / boundary, half_width, half_height);
}

Cairo::Matrix GraphArea::canvas_to_graph() const {
    Gtk::Allocation allocation = get_allocation();

    const double half_width = 0.5 * allocation.get_width(), half_height = 0.5 * allocation.get_height();
    const double boundary = gui_.creativity_->boundary();
    // The second value below is negative because higher y canvas values are down the screen, which
    // are *lower* graph values.
    //
    // After scaling, the (x+, y+) screen location will be scaled into ([0,2*b], [0,-2*b])
    // coordinates, which we need to shift into ([-b,b], [-b,b]) by subtracting b from x and adding
    // b to y.
    return Cairo::Matrix(boundary / half_width, 0, 0, -boundary / half_height, -boundary, boundary);
}

bool GraphArea::on_draw(const Cairo::RefPtr<Cairo::Context> &cr_grapharea) {
    Gtk::Allocation allocation = get_allocation();
    const int width = allocation.get_width();
    const int height = allocation.get_height();
    const double boundary = gui_.creativity_->boundary();

    if (gui_.state_num_ == 0) {
        // No states at all: just draw a blank screen
        cr_grapharea->set_source(design.colour.background);
        cr_grapharea->paint();
        return true;
    }

    if (drawing_cache_width_ != width || drawing_cache_height_ != height) {
        // If the width or height has changed from the cached image, clear the cache.
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
        cr->set_source(design.colour.background);
        cr->paint();

        auto trans = graph_to_canvas();

        // Things get draw here in order from least to most important (since later things are drawn
        // on top of earlier things).  This means we have to loop over the same vector a few times,
        // unfortunately.

        // Draw lines from readers to newly purchased books (unless we're not showing newly
        // purchased books)
        if (design.enabled.reading) {
            cr->save();
            cr->set_source(design.colour.reading);
            cr->set_dash(design.dash.reading, 0);
            cr->set_line_width(design.stroke_width.reading);
            for (auto &rpair : state->readers) {
                auto &r = rpair.second;

                for (auto &book_id : r.new_books) {
                    auto &b = state->books.at(book_id);
                    drawWrappingLine(cr, trans, r.position, b.position);
                }
            }
            cr->stroke();
            cr->restore();
        }

        // Draw friendship/sharing links
        if (design.enabled.friendship) {
            cr->save();
            cr->set_source(design.colour.friendship);
            cr->set_dash(design.dash.friendship, 0);
            cr->set_line_width(design.stroke_width.friendship);
            for (auto &rpair : state->readers) {
                auto &r = rpair.second;
                for (auto &fid : r.friends) {
                    if (r.id < fid) { // lines are symmetric: only draw each link once
                        auto f = state->readers.find(fid);
                        if (f != state->readers.end())
                            drawWrappingLine(cr, trans, r.position, f->second.position);
                    }
                }
            }
            cr->stroke();
            cr->restore();
        }


        // Draw books.  On-market and off-market books have different colours, and newer books are
        // larger than older books.
        if (design.enabled.book_live or design.enabled.book_dead) {
            for (auto &bpair : state->books) {
                auto &b = bpair.second;
                // Give new books a larger cross: a brand new book gets a point scaled larges; this
                // scaling decreases linearly until reading the regular size.

                if (b.market) {
                    if (design.enabled.book_live) {
                        const double scale = std::max(1.0, design.size.book_live_scale_a - design.size.book_live_scale_b*b.age);
                        drawPoint(cr, trans, b.position[0], b.position[1], design.style.book_live,
                                design.colour.book_live, design.size.book_live*scale, design.stroke_width.book_live);
                    }
                }
                else {
                    if (design.enabled.book_dead) {
                        const double scale = std::max(1.0, design.size.book_dead_scale_a - design.size.book_dead_scale_b*b.age);
                        drawPoint(cr, trans, b.position[0], b.position[1], design.style.book_dead,
                                design.colour.book_dead, design.size.book_dead*scale, design.stroke_width.book_dead);
                    }
                }
            }
        }

        // Readers have utility circles, indicating how much about the base 1000 utility the reader
        // was in the period.  (Don't draw anything in the initialization period, though, since
        // everyone just has 0 utility then).
        if ((design.enabled.utility_gain or design.enabled.utility_loss) and state->t > 0) {
            cr->save();
            cr->set_dash(design.dash.utility, 0);
            cr->set_line_width(design.stroke_width.utility);
            for (auto &rpair : state->readers) {
                auto &r = rpair.second;
                double rx = r.position[0], ry = r.position[1];
                double radius = 0.0;
                Colour colour = design.colour.utility_gain;
                if (r.u > 1000) { /// 1000 is the starting utility
                    radius = design.size.utility_gain_scale * std::log(r.u - 999);
                }
                else if (r.u < 1000) {
                    radius = design.size.utility_loss_scale * std::log(1001 - r.u);
                    colour = design.colour.utility_loss;
                }
                ERIS_DBGVAR(radius);

                if (radius > 0.0)
                    drawCanvasCircle(cr, trans, rx, ry, radius, colour, design.stroke_width.utility, design.stroke_width.utility_radial);
            }
            cr->restore();
        }

        // Lines from each book to its author
        if (design.enabled.author_live or design.enabled.author_dead) {
            cr->save();
            for (auto &bpair : state->books) {
                auto &b = bpair.second;
                if (b.market) {
                    if (not design.enabled.author_live) continue;

                    cr->set_source(design.colour.author_live);
                    cr->set_dash(design.dash.author_live, 0);
                    cr->set_line_width(design.stroke_width.author_live);
                }
                else {
                    if (not design.enabled.author_dead) continue;

                    cr->set_source(design.colour.author_dead);
                    cr->set_dash(design.dash.author_dead, 0);
                    cr->set_line_width(design.stroke_width.author_dead);
                }
                drawWrappingLine(cr, trans, b.position, state->readers.at(b.author).position);
                cr->stroke();
            }
            cr->restore();
        }

        // Draw readers
        if (design.enabled.reader) {
            for (auto &rpair : state->readers) {
                auto &r = rpair.second;
                double rx = r.position[0], ry = r.position[1];
                // Draw the reader
                drawPoint(cr, trans, rx, ry,
                        design.style.reader, design.colour.reader, design.size.reader, design.stroke_width.reader);
            }
        }

        // Add axes (last, so that they are on top)
        double zero_x = 0, zero_y = 0;
        trans.transform_point(zero_x, zero_y);

        // Round the transformed location to the nearest value; it'll be ever so slighly inaccurate,
        // but will give a double pixel line instead of a slightly blurry line (from partial pixel coverage).
        zero_x = std::round(zero_x);
        zero_y = std::round(zero_y);

        cr->set_line_width(design.stroke_width.axes);
        cr->set_source(design.colour.axes);
        cr->move_to(zero_x, 0);
        cr->line_to(zero_x, height);
        cr->stroke();
        cr->move_to(0, zero_y);
        cr->line_to(width, zero_y);
        cr->stroke();

        // Tick marks
        int tick_num = 0;
        for (double tickpos = design.style.tick_every; tickpos <= boundary; tickpos += design.style.tick_every) {
            const bool big_tick = ++tick_num % design.style.tick_big == 0;
            const double curr_tick_size = big_tick ? design.size.tick_big : design.size.tick;
            cr->set_line_width(big_tick ? design.stroke_width.axes_ticks_big : design.stroke_width.axes_ticks);

            // x-axis ticks:
            for (double x : {tickpos,-tickpos}) {
                double y = 0;
                trans.transform_point(x, y);
                x = std::round(x + 0.5) - 0.5; // Round to nearest 0.5 for crisper lines
                cr->move_to(x, y - 0.5*curr_tick_size);
                cr->rel_line_to(0, curr_tick_size);
            }
            // y-axis ticks:
            for (double y : {tickpos,-tickpos}) {
                double x = 0;
                trans.transform_point(x, y);
                y = std::round(y + 0.5) - 0.5; // Round to nearest 0.5 for crisper lines
                cr->move_to(x - 0.5*curr_tick_size, y);
                cr->rel_line_to(curr_tick_size, 0);
            }

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

void GraphArea::drawPoint(
        const Cairo::RefPtr<Cairo::Context> &cr, const Cairo::Matrix &trans, double x, double y,
        const PointType &type, const Colour &colour, const double radius, const double width, bool virt) {

    const double boundary = gui_.creativity_->boundary();

    // Before we transform x and y, consider whether there are wrapped versions of the point we need
    // to worry about; if so, draw them first.
    if (not virt) {
        double dim_x = 2*boundary,
               dim_y = 2*boundary;
        bool wrap_right = x + radius > boundary,
             wrap_left  = x - radius < -boundary,
             wrap_up    = y + radius > boundary,
             wrap_down  = y - radius < -boundary;

        // There are 8 wrapping possibilities to consider; the first 4 are corners, the last 4 are
        // single-dimension edge wraps.  Corners are just like a wrap in each wrapped edge, but
        // *also* have to wrap to the opposite corner.
        //
        // Note also that this doesn't properly handle cases where radius > 2*boundary: such would
        // require wrapping multiple times, or on opposite sides at once (i.e. right *and* left
        // wraps).  Since that shouldn't happen in practice, this isn't a big limitation.
        if (wrap_left and wrap_up) // Upper left
            drawPoint(cr, trans, x+dim_x, y-dim_y, type, colour, radius, width, true); // Lower right
        else if (wrap_left and wrap_down) // Lower left
            drawPoint(cr, trans, x+dim_x, y+dim_y, type, colour, radius, width, true); // Upper right
        else if (wrap_right and wrap_up) // Upper right
            drawPoint(cr, trans, x-dim_x, y-dim_y, type, colour, radius, width, true); // Lower left
        else if (wrap_right and wrap_down) // Lower right
            drawPoint(cr, trans, x-dim_x, y+dim_y, type, colour, radius, width, true); // Upper left

        // Now the edge mirroring
        if (wrap_left)
            drawPoint(cr, trans, x+dim_x, y, type, colour, radius, width, true);
        else if (wrap_right)
            drawPoint(cr, trans, x-dim_x, y, type, colour, radius, width, true);

        if (wrap_up)
            drawPoint(cr, trans, x, y-dim_y, type, colour, radius, width, true);
        else if (wrap_down)
            drawPoint(cr, trans, x, y+dim_y, type, colour, radius, width, true);
    }

    // Get the user-space coordinates of the point
    trans.transform_point(x, y);

    cr->save();
    cr->set_line_width(width);
    cr->set_source(colour);

    switch (type) {
        case PointType::X:
            {
                // We're drawing 45-degree lines, so get the linear distance on each dimension:
                double edge = radius * sqrt(2.0);

                cr->move_to(x-0.5*edge, y-0.5*edge); // Upper left
                cr->rel_line_to(edge, edge); // ... line to bottom right

                cr->rel_move_to(-edge, 0); // Lower left
                cr->rel_line_to(edge, -edge); // ... to upper right
                break;
            }
        case PointType::SQUARE:
            {
                // The SQUARE vertices are the the same as the X type, above
                double edge = radius * sqrt(2.0);
                cr->rectangle(x-0.5*edge, y-0.5*edge, edge, edge);
                break;
            }
        case PointType::CROSS:
            {
                // Easier than the above: the lines are simply of length 2*radius
                // Horizontal line:
                cr->move_to(x - radius, y);
                cr->rel_line_to(2*radius, 0);
                // Verical line:
                cr->move_to(x, y - radius);
                cr->rel_line_to(0, 2*radius);
                break;
            }
    }

    cr->stroke();
    cr->restore();
}

void GraphArea::drawGraphCircle(const Cairo::RefPtr<Cairo::Context> &cr, const Cairo::Matrix &trans,
        double cx, double cy, double r, const Colour &colour, double stroke_width, double radial_stroke_width) {
    // Apply the full transformation: then we can just draw the circle in *graph* space which,
    // post-transformation, will be exactly the desired oval.
    cr->save();
    cr->transform(trans);

    // Set up the circle in the transformed (graph) coordinates, but *don't* draw it yet (if we do,
    // the line width will also be transformed, which we don't want--it would be too thin on the
    // more compressed axis).
    cr->arc(cx, cy, r, 0.0, 2*M_PI);

    // Undo the transformation, then change the color and line width and actually draw the circle
    cr->restore();
    cr->save();
    cr->set_line_width(stroke_width);
    cr->set_source(colour);

    cr->stroke();
    cr->restore();

    if (radial_stroke_width > 0) {
        cr->save();
        cr->transform(trans);
        // Draw a radial line at a quasi-random angle (based on the position and radius values) to the
        // center of the circle.  The quasi-random calculation means the same circle gets the same angle
        // when redrawn.
        double angle = std::fmod(1000000 * cx * cy * r, 2*M_PI);
        cr->move_to(cx, cy);
        cr->rel_line_to(r*std::sin(angle), r*std::cos(angle));
        cr->restore();
        cr->save();
        cr->set_line_width(radial_stroke_width);
        cr->set_source(colour);
        cr->stroke();
        cr->restore();
    }
}

void GraphArea::drawCanvasCircle(const Cairo::RefPtr<Cairo::Context> &cr, const Cairo::Matrix &trans,
        double cx, double cy, double r, const Colour &colour, double stroke_width, double radial_stroke_width) {

    // Draw a radial line at a quasi-random angle (based on the position and radius values) to the
    // center of the circle.  The quasi-random calculation means the same circle gets the same angle
    // when redrawn.  Calculate this *before* transforming the point so that the angle remains the
    // same across canvas resizing.
    const double radial_angle = std::fmod(1000000 * cx * cy * r, 2*M_PI);

    trans.transform_point(cx, cy);

    cr->save();
    cr->arc(cx, cy, r, 0.0, 2*M_PI);
    cr->set_line_width(stroke_width);
    cr->set_source(colour);
    cr->stroke();

    // Now draw the radial line (if non-zero):
    if (radial_stroke_width > 0) {
        cr->set_line_width(radial_stroke_width);

        cr->move_to(cx, cy);
        cr->rel_line_to(r*std::sin(radial_angle), r*std::cos(radial_angle));
        cr->set_line_width(radial_stroke_width);
        cr->stroke();
    }

    cr->restore();
}

void GraphArea::resetCache() {
    drawing_cache_width_ = -1;
}

} }

