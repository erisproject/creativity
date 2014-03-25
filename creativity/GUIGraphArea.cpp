#include "creativity/GUIGraphArea.hpp"
#include "creativity/GUI.hpp"
#include "creativity/Reader.hpp"
#include "creativity/Book.hpp"

namespace creativity {

GUIGraphArea::GUIGraphArea(const double &top, const double &right, const double &bottom, const double &left,
        GUI &gui, std::shared_ptr<eris::Simulation> sim) :
    gui_{gui}, sim_{sim}, bounds_{{top,right,bottom,left}}
{}

bool GUIGraphArea::on_draw(const Cairo::RefPtr<Cairo::Context> &cr) {

    Gtk::Allocation allocation = get_allocation();
    const int width = allocation.get_width();
    const int height = allocation.get_height();

    // The second of these is typically negative, which is correct because "up" in the graph
    // translates to a lower coordinate on the screen (since positive y coordinates are down the
    // screen).
    const double gwidth = bounds_[RIGHT] - bounds_[LEFT],
          gheight = bounds_[BOTTOM] - bounds_[TOP];

    cr->save();

    // Lock out the simulation until we finish redrawing
    auto lock = sim_->runLock();

    // Paint the background white
    cr->set_source_rgb(1,1,1);
    cr->paint();

    // Build a transformation that converts from positions to canvas coordinates
    auto trans = Cairo::identity_matrix();
    trans.scale(width / gwidth, height / gheight);
    trans.translate(-bounds_[LEFT], -bounds_[TOP]);

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
    const double tick_max_x = std::max(std::abs(bounds_[RIGHT]), std::abs(bounds_[LEFT]));
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
    const double tick_max_y = std::max(std::abs(bounds_[TOP]), std::abs(bounds_[BOTTOM]));
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

    cr->restore();

    for (auto &r : sim_->agents<Reader>()) {
        drawPoint(cr, trans, r->position()[0], r->position()[1], PointType::X);
    }
    for (auto &b : sim_->goods<Book>()) {
        // Change the radius of the book: new books have radius 1 (in graph space), the radius decreases each period
        // by 0.05 until 0, at which point we stop drawing the circle.
        auto radius = 1.0 - b->age() / 20.0;

        drawPoint(cr, trans, b->position()[0], b->position()[1], PointType::SQUARE);
        if (radius > 0) {
            drawCircle(cr, trans, b->position()[0], b->position()[1], radius, CircleType::A);
        }
    }
    return true;
}

void GUIGraphArea::drawPoint(const Cairo::RefPtr<Cairo::Context> &cr, const Cairo::Matrix &trans, double x, double y, const PointType &type) {
    // Get the user-space coordinates of the point
    trans.transform_point(x, y);

    cr->save();
    cr->set_line_width(2.0);

    switch (type) {
        case PointType::X:
            {
                // We're drawing 45-degree lines, so get the linear distance on each dimension:
                double edge = gui_.pointMarkerRadius * sqrt(2.0);

                cr->set_source_rgb(1.0, 0, 0); // red

                cr->move_to(x-0.5*edge, y-0.5*edge); // Upper left
                cr->rel_line_to(edge, edge); // ... line to bottom right

                cr->rel_move_to(-edge, 0); // Lower left
                cr->rel_line_to(edge, -edge); // ... to upper right
                break;
            }
        case PointType::SQUARE:
            {
                // The SQUARE vertices are the the same as the X type, above
                double edge = gui_.pointMarkerRadius * sqrt(2.0);

                cr->set_source_rgb(0.5, 0, 1.0); // Purple
                cr->rectangle(x-0.5*edge, y-0.5*edge, edge, edge);
                break;
            }
        case PointType::CROSS:
            {
                cr->set_source_rgb(0.0, 0.25, 1.0); // Blue with a hint of green
                // Horizontal line:
                cr->move_to(x - gui_.pointMarkerRadius, y);
                cr->rel_line_to(2*gui_.pointMarkerRadius, 0);

                // Verical line:
                cr->move_to(x, y - gui_.pointMarkerRadius);
                cr->rel_line_to(0, 2*gui_.pointMarkerRadius);
                break;
            }
    }

    cr->stroke();
    cr->restore();
}

void GUIGraphArea::drawCircle(const Cairo::RefPtr<Cairo::Context> &cr, const Cairo::Matrix &trans,
        const double &cx, const double &cy, const double &r, const CircleType &type) {
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
    cr->set_line_width(2.0);
    switch (type) {
        case CircleType::A:
            cr->set_source_rgba(1.0, 0.55, 0, 0.5); // Orange
            break;
        case CircleType::B:
            cr->set_source_rgba(0.133, 0.545, 0.133, 0.5);
            break;
    }

    cr->stroke();
    cr->restore();
}

}

