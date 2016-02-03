#include "creativity/gui/GraphArea.hpp"
#include "creativity/gui/GUI.hpp"
#include "creativity/state/Storage.hpp"
#include "creativity/state/State.hpp"
#include "creativity/Creativity.hpp"
#include <boost/math/constants/constants.hpp>
#include <cairomm/context.h>
#include <cairomm/matrix.h>
#include <cairomm/pattern.h>
#include <cairomm/refptr.h>
#include <cairomm/surface.h>
#include <cairomm/enums.h>
#include <gtkmm/widget.h>
#include <cmath>
#include <algorithm>
#include <map>
#include <utility>
#include <memory>


using namespace eris;
using namespace creativity::state;

namespace creativity { namespace gui {

GraphArea::GraphArea(GUI &gui) :
    gui_(gui), wpb_({0.0,0.0})
{}

Cairo::Matrix GraphArea::graph_to_canvas() const {
    Gtk::Allocation allocation = get_allocation();

    const double half_width = 0.5 * allocation.get_width(), half_height = 0.5 * allocation.get_height();
    const double boundary = gui_.creativity_.parameters.boundary;
    // The second of these is negative because higher y canvas values are down the screen, which are
    // *lower* graph values.
    return Cairo::Matrix(half_width / boundary, 0, 0, -half_height / boundary, half_width, half_height);
}

Cairo::Matrix GraphArea::canvas_to_graph() const {
    Gtk::Allocation allocation = get_allocation();

    const double half_width = 0.5 * allocation.get_width(), half_height = 0.5 * allocation.get_height();
    const double boundary = gui_.creativity_.parameters.boundary;
    // The second value below is negative because higher y canvas values are down the screen, which
    // are *lower* graph values.
    //
    // After scaling, the (x+, y+) screen location will be scaled into ([0,2*b], [0,-2*b])
    // coordinates, which we need to shift into ([-b,b], [-b,b]) by subtracting b from x and adding
    // b to y.
    return Cairo::Matrix(boundary / half_width, 0, 0, -boundary / half_height, -boundary, boundary);
}

Position GraphArea::graph_position(const Position &p) {
    if (p.dimensions == 2) return p;
    if (p.dimensions == 1) {
        const double r = 0.9 * gui_.creativity_.parameters.boundary;
        const double angle = p[0] / gui_.creativity_.parameters.boundary * boost::math::constants::pi<double>();
        double x = r * std::cos(angle);
        double y = r * std::sin(angle);
        return Position({x, y});
    }
    // 3+ dimensions
    return p.subdimensions({design.subdimensionX, design.subdimensionY});
}

bool GraphArea::on_draw(const Cairo::RefPtr<Cairo::Context> &cr_grapharea) {
    Gtk::Allocation allocation = get_allocation();
    const int width = allocation.get_width();
    const int height = allocation.get_height();

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

        std::shared_ptr<const State> state, prev_state;

        {
            auto st = gui_.creativity_.storage();
            if (st.first->empty())
                state = std::make_shared<State>();
            else {
                state = (*st.first)[gui_.state_curr_];
                if (gui_.state_curr_ > 0) prev_state = (*st.first)[gui_.state_curr_ - 1];
            }
        }

        cr->save();

        // Paint the background white
        cr->set_source(design.colour.background);
        cr->paint();

        auto trans = graph_to_canvas();


        // In 1D mode, draw the circle first (in 2+D, we draw the axes at the end so that they're on
        // top of everything).
        if (gui_.creativity_.parameters.dimensions == 1 and design.enabled.axes) {
            drawCircularAxis(cr, trans);
        }

        // Things get draw here in order from least to most important (since later things are drawn
        // on top of earlier things).  This means we have to loop over the same vector a few times,
        // unfortunately.

        // Draw reader movement lines (if the reader existed in the previous period).
        if (design.enabled.reader and design.enabled.movement and prev_state) {
            cr->save();
            cr->set_dash(design.dash.movement, 0);
            cr->set_line_width(design.stroke_width.movement);
            double r, g, b, a;
            design.colour.movement->get_rgba(r,g,b,a);
            a = std::max(0.0, std::min(1.0, a * design.style.movement_alpha_multiplier));
            auto start_colour = Cairo::SolidPattern::create_rgba(r, g, b, a);

            for (auto &rpair : state->readers) {
                if (prev_state->readers.count(rpair.first) > 0) {
                    auto &was_at = prev_state->readers.at(rpair.first).position;
                    drawLine(cr, trans, was_at, rpair.second.position, 5.0, start_colour, design.colour.movement);
                }
            }
            cr->restore();
        }

        // Draw lines from readers to newly purchased books (unless we're not showing newly
        // purchased books)
        if (design.enabled.reading and design.enabled.reader and
                (design.enabled.book_live or design.enabled.book_dead or design.enabled.book_public)) {
            cr->save();
            cr->set_source(design.colour.reading);
            cr->set_dash(design.dash.reading, 0);
            cr->set_line_width(design.stroke_width.reading);
            for (auto &rpair : state->readers) {
                auto &r = rpair.second;

                for (auto &book_id : r.new_books) {
                    auto &b = state->books.at(book_id);
                    bool show = (b.market_private ? design.enabled.book_live :
                            b.market_public() ? design.enabled.book_public :
                            design.enabled.book_dead);
                    if (show)
                        drawLine(cr, trans, r.position, b.position);
                }
            }
            cr->stroke();
            cr->restore();
        }

        // Draw friendship/sharing links
        if (design.enabled.friendship and design.enabled.reader) {
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
                            drawLine(cr, trans, r.position, f->second.position);
                    }
                }
            }
            cr->stroke();
            cr->restore();
        }


        // Draw books.  On-market and off-market books have different colours, and newer books are
        // larger than older books.
        if (design.enabled.book_live or design.enabled.book_dead or design.enabled.book_public) {
            for (auto &bpair : state->books) {
                auto &b = bpair.second;
                // Give new books a larger cross: a brand new book gets a point scaled larges; this
                // scaling decreases linearly until reading the regular size.
                const double size = std::max(1.0, design.size.book_scale_a - design.size.book_scale_b*(state->t - b.created))
                    * design.size.book;

                if (b.market_private) {
                    if (design.enabled.book_live) {
                        drawPoint(cr, trans, b.position, design.style.book_live,
                                design.colour.book_live, size, design.stroke_width.book_live);
                    }
                }
                else if (b.market_public()) {
                    if (design.enabled.book_public) {
                        drawPoint(cr, trans, b.position, design.style.book_public,
                                design.colour.book_public, size, design.stroke_width.book_public);
                    }
                }
                else {
                    if (design.enabled.book_dead) {
                        drawPoint(cr, trans, b.position, design.style.book_dead,
                                design.colour.book_dead, size, design.stroke_width.book_dead);
                    }
                }
            }
        }

        // Readers have utility circles, indicating how much about the base 1000 utility the reader
        // was in the period.  (Don't draw anything in the initialization period, though, since
        // everyone just has 0 utility then).
        if ((design.enabled.utility_gain or design.enabled.utility_loss) and design.enabled.reader and state->t > 0) {
            cr->save();
            cr->set_dash(design.dash.utility, 0);
            cr->set_line_width(design.stroke_width.utility);
            for (auto &rpair : state->readers) {
                auto &r = rpair.second;
                double radius = 0.0;
                Colour colour = design.colour.utility_gain;
                if (r.u > 1000) { /// 1000 is the starting utility
                    radius = design.size.utility_gain_scale * std::log(r.u - 999);
                }
                else if (r.u < 1000) {
                    radius = design.size.utility_loss_scale * std::log(1001 - r.u);
                    colour = design.colour.utility_loss;
                }

                if (radius > 0.0)
                    drawCanvasCircle(cr, trans, r.position, radius, colour, design.stroke_width.utility, design.stroke_width.utility_radial);
            }
            cr->restore();
        }

        // Lines from each book to its author
        bool show_author_on = design.enabled.reader and design.enabled.author_live and design.enabled.book_live,
             show_author_off = design.enabled.reader and design.enabled.author_dead and design.enabled.book_dead,
             show_author_pub = design.enabled.reader and design.enabled.author_public and design.enabled.book_public;

        if (show_author_on or show_author_off or show_author_pub) {
            cr->save();
            for (auto &bpair : state->books) {
                auto &b = bpair.second;
                if (b.market_private) {
                    if (not show_author_on) continue;

                    cr->set_source(design.colour.author_live);
                    cr->set_dash(design.dash.author_live, 0);
                    cr->set_line_width(design.stroke_width.author_live);
                }
                else if (b.market_public()) {
                    if (not show_author_pub) continue;

                    cr->set_source(design.colour.author_public);
                    cr->set_dash(design.dash.author_public, 0);
                    cr->set_line_width(design.stroke_width.author_public);
                }
                else {
                    if (not show_author_off) continue;

                    cr->set_source(design.colour.author_dead);
                    cr->set_dash(design.dash.author_dead, 0);
                    cr->set_line_width(design.stroke_width.author_dead);
                }
                drawLine(cr, trans, b.position, state->readers.at(b.author).position);
                cr->stroke();
            }
            cr->restore();
        }

        // Draw readers
        if (design.enabled.reader) {
            for (auto &rpair : state->readers) {
                auto &r = rpair.second;
                // Draw the reader
                drawPoint(cr, trans, r.position,
                        design.style.reader, design.colour.reader, design.size.reader, design.stroke_width.reader);
            }
        }

        // Add axes (last, so that they are on top)
        if (gui_.creativity_.parameters.dimensions != 1 and design.enabled.axes) {
            drawAxes(cr, trans);
        }
    }

    cr_grapharea->set_source(drawing_cache_[gui_.state_curr_], 0, 0);
    cr_grapharea->paint();

    return true;
}

void GraphArea::drawLine(const Cairo::RefPtr<Cairo::Context> &cr, const Cairo::Matrix &trans, const Position &from_raw, const Position &to_raw,
        double min_length, Colour start_colour, Colour end_colour) {

    Position from = graph_position(from_raw), to = graph_position(to_raw);

    Position v; // The line vector
    const double &boundary = gui_.creativity_.parameters.boundary;
    const double bshift = 2*boundary;

    // If we're doing a 1D graph, this is much simpler: we don't have to worry about extra lines,
    // wrapping edges, etc.: all lines are direct (across the circle).
    if (from_raw.dimensions == 1) {
        v = to - from;
    }
    else {
        // Otherwise we have to wrap
        if (not wpb_.wrapped(0)) {
            wpb_ = WrappedPositionalBase({0.0,0.0}, {-boundary, -boundary}, {boundary, boundary});
        }
        wpb_.moveTo(from);
        // Get the shortest-path vector to the target, first in graph vectors, then translated to canvas
        // vectors:
        v = wpb_.vectorTo(to);
    }
    // Translate the vector into canvas space:
    double v_canvas[2] = {v[0], v[1]};
    trans.transform_distance(v_canvas[0], v_canvas[1]);

    if (std::hypot(v_canvas[0], v_canvas[1]) < min_length)
        return;

    // Build a list of all the "from" points we actually need to draw from.  Obviously we need a
    // vector from the original one, but we may also need to draw a vector from 1 or 3 virtual
    // source locations (i.e. locations in wrapped dimensions).
    std::vector<std::pair<double,double>> source_locations;
    source_locations.emplace_back(from[0], from[1]);

    if (from_raw.dimensions > 1) {
        // Figure out if we're crossing any boundaries:
        bool wrap_right = (v[0] > 0 and to[0] < from[0]),
             wrap_left = (v[0] < 0 and to[0] > from[0]),
             wrap_top = (v[1] > 0 and to[1] < from[1]),
             wrap_bottom = (v[1] < 0 and to[1] > from[1]);

        if (wrap_right) {
            source_locations.emplace_back(from[0] - bshift, from[1]); // Left virtual
            // Check for corner virtuals:
            if (wrap_top)
                source_locations.emplace_back(from[0] - bshift, from[1] - bshift); // Bottom-left virtual
            else if (wrap_bottom)
                source_locations.emplace_back(from[0] - bshift, from[1] + bshift); // Top-left virtual
        }
        else if (wrap_left) {
            source_locations.emplace_back(from[0] + bshift, from[1]); // Right virtual
            // Check for corner virtuals:
            if (wrap_top)
                source_locations.emplace_back(from[0] + bshift, from[1] - bshift); // Bottom-right virtual
            else if (wrap_bottom)
                source_locations.emplace_back(from[0] + bshift, from[1] + bshift); // Top-right virtual
        }
        if (wrap_top)
            source_locations.emplace_back(from[0], from[1] - bshift); // Bottom virtual
        else if (wrap_bottom)
            source_locations.emplace_back(from[0], from[1] + bshift); // Top virtual
    }

    // Now figure out the gradient/colour stuff, if requested
    cr->save();
    struct rgba { double r, g, b, a; };
    struct { bool active = false; rgba start; rgba end; } gradient;
    if (start_colour) {
        if (end_colour) {
            gradient.active = true;
            start_colour->get_rgba(gradient.start.r, gradient.start.g, gradient.start.b, gradient.start.a);
            end_colour->get_rgba(gradient.end.r, gradient.end.g, gradient.end.b, gradient.end.a);
        }
        else {
            cr->set_source(start_colour);
        }
    }

    // else start_colour is unset, so don't touch the colour
    // source_locations now holds all the virtual line segment starting points (in graph notation),
    // from which we need to draw the directional vector `v_canvas` (in canvas notation).
    for (auto &s : source_locations) {
        trans.transform_point(s.first, s.second); // Convert to canvas coordinate
        if (gradient.active) {
            auto grad = Cairo::LinearGradient::create(s.first, s.second, s.first + v_canvas[0], s.second + v_canvas[1]);
            grad->add_color_stop_rgba(0.0, gradient.start.r, gradient.start.g, gradient.start.b, gradient.start.a);
            grad->add_color_stop_rgba(1.0, gradient.end.r, gradient.end.g, gradient.end.b, gradient.end.a);
            cr->set_source(grad);
        }
        cr->move_to(s.first, s.second);
        cr->rel_line_to(v_canvas[0], v_canvas[1]);
        cr->stroke();
    }

    cr->restore(); // Undo any colour changes
}

void GraphArea::drawPoint(
        const Cairo::RefPtr<Cairo::Context> &cr, const Cairo::Matrix &trans, const Position &point,
        const PointType &type, const Colour &colour, const double radius, const double width) {

    const double &boundary = gui_.creativity_.parameters.boundary;

    Position graphpos = graph_position(point);
    const double &x = graphpos[0];
    const double &y = graphpos[1];

    // Before we transform x and y, consider whether there are wrapped versions of the point we need
    // to worry about; if so, draw them first.
    const double bshift = 2*boundary;
    bool wrap_right = x + radius > boundary,
         wrap_left  = x - radius < -boundary,
         wrap_up    = y + radius > boundary,
         wrap_down  = y - radius < -boundary;

    // There are 8 wrapping possibilities to consider; the first 4 are corners, the last 4 are
    // single-dimension edge wraps.  Corners are just like a wrap in each wrapped edge, but *also*
    // have to wrap to the opposite corner.
    //
    // Note also that this doesn't properly handle cases where radius > 2*boundary: such would
    // require wrapping multiple times, or on opposite sides at once (i.e. right *and* left wraps).
    // Since that shouldn't happen in practice, this isn't a big limitation.
    if (wrap_left and wrap_up) // Upper left
        drawPointSingle(cr, trans, x+bshift, y-bshift, type, colour, radius, width); // Lower right
    else if (wrap_left and wrap_down) // Lower left
        drawPointSingle(cr, trans, x+bshift, y+bshift, type, colour, radius, width); // Upper right
    else if (wrap_right and wrap_up) // Upper right
        drawPointSingle(cr, trans, x-bshift, y-bshift, type, colour, radius, width); // Lower left
    else if (wrap_right and wrap_down) // Lower right
        drawPointSingle(cr, trans, x-bshift, y+bshift, type, colour, radius, width); // Upper left

    // Now the edge mirroring
    if (wrap_left)
        drawPointSingle(cr, trans, x+bshift, y, type, colour, radius, width);
    else if (wrap_right)
        drawPointSingle(cr, trans, x-bshift, y, type, colour, radius, width);

    if (wrap_up)
        drawPointSingle(cr, trans, x, y-bshift, type, colour, radius, width);
    else if (wrap_down)
        drawPointSingle(cr, trans, x, y+bshift, type, colour, radius, width);

    drawPointSingle(cr, trans, x, y, type, colour, radius, width);
}


void GraphArea::drawPointSingle(
        const Cairo::RefPtr<Cairo::Context> &cr, const Cairo::Matrix &trans, double x, double y,
        const PointType &type, const Colour &colour, const double radius, const double width) {

    // Get the canvas coordinates of the point
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
        const eris::Position &centre, double r, const Colour &colour, double stroke_width, double radial_stroke_width) {

    Position c = graph_position(centre);
    double cx = c[0], cy = c[1];

    // Draw a radial line at a quasi-random angle (based on the position and radius values) to the
    // center of the circle.  The quasi-random calculation means the same circle gets the same angle
    // when redrawn.  Also calculate this *before* transforming the point so that the angle remains
    // the same across canvas resizing.
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

void GraphArea::drawCircularAxis(const Cairo::RefPtr<Cairo::Context> &cr, const Cairo::Matrix &trans) {
    double zero_x = 0, zero_y = 0;
    trans.transform_point(zero_x, zero_y);

    const double &boundary = gui_.creativity_.parameters.boundary;
    const double r = 0.9 * boundary;
    drawGraphCircle(cr, trans, 0, 0, r, design.colour.axes, design.stroke_width.axes);

    cr->save();
    cr->set_line_width(design.stroke_width.axes);
    cr->set_source(design.colour.axes);

    // Tick marks
    int tick_num = 0;
    for (double tick = 0; tick <= boundary; tick += design.style.tick_every) {
        const bool big_tick = tick_num++ % design.style.tick_big == 0;
        const double curr_tick_size = big_tick ? design.size.tick_big : design.size.tick;
        cr->set_line_width(big_tick ? design.stroke_width.axes_ticks_big : design.stroke_width.axes_ticks);
        for (const auto &sign : {1, -1}) {
            if (sign < 0 and tick == 0) continue;

            // Get the graph location of the tick we're interested in
            Position tickpos = graph_position({sign*tick});
            // using tickpos directly yields the vector perpendicular (in graph space) to the circle
            // at the location of the tick.  We want to transform the perpendicular of *that* (and
            // so get the tangent to the circle):
            double tangdx = tickpos[1], tangdy = -tickpos[0];
            // Now transform the tangent into canvas space:
            trans.transform_distance(tangdx, tangdy);
            // So we now have a tangent of the transformed circle on the canvas; get the
            // perpendicular of that for our tick direction
            double tickx = tangdy, ticky = -tangdx;
            // Now figure out how much to scale it to get a vector of the desired tick length
            double scale = curr_tick_size / std::hypot(tickx, ticky);
            tickx *= scale;
            ticky *= scale;

            // Get the tickpos translated into canvas coordinates:
            double tx = tickpos[0], ty = tickpos[1];
            trans.transform_point(tx, ty);
            // Back up half a tick from there, then draw forward a full tick:
            cr->move_to(tx - 0.5*tickx, ty - 0.5*ticky);
            cr->rel_line_to(tickx, ticky);
        }
        cr->stroke();
    }

    cr->restore();

}

void GraphArea::drawAxes(const Cairo::RefPtr<Cairo::Context> &cr, const Cairo::Matrix &trans) {
    double zero_x = 0, zero_y = 0;
    trans.transform_point(zero_x, zero_y);

    // Round the transformed location to the nearest value; it'll be ever so slighly inaccurate,
    // but will give a double pixel line instead of a slightly blurry line (from partial pixel coverage).
    zero_x = std::round(zero_x);
    zero_y = std::round(zero_y);

    cr->save();

    cr->set_line_width(design.stroke_width.axes);
    cr->set_source(design.colour.axes);
    cr->move_to(zero_x, 0);
    cr->line_to(zero_x, drawing_cache_height_);
    cr->stroke();
    cr->move_to(0, zero_y);
    cr->line_to(drawing_cache_width_, zero_y);
    cr->stroke();

    const double &boundary = gui_.creativity_.parameters.boundary;

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

void GraphArea::resetCache() {
    drawing_cache_width_ = -1;
}

} }

