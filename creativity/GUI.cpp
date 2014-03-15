#include "creativity/GUI.hpp"
#include "creativity/config.hpp"
#include <iostream>
#include <cstdlib>

namespace creativity {

template <class T>
std::unique_ptr<T> builder_widget(std::string widget_name, Glib::RefPtr<Gtk::Builder> builder) {
    T *widget;
    builder->get_widget(widget_name, widget);
    return std::unique_ptr<T>(widget);
}

GUI::Exception::Exception(const std::string &what) : std::runtime_error(what) {}

GUI::GUI() {}

void GUI::start(int argc, char *argv[]) {
    auto app = Gtk::Application::create(argc, argv, "ca.imaginary.test.cairo-drawing",
            Gio::APPLICATION_NON_UNIQUE);
    auto builder = Gtk::Builder::create();

    std::string datadir(DATADIR);
    char *envdatadir = getenv("CREATIVITY_DATADIR");
    if (envdatadir) {
        std::string stddatadir(envdatadir);
        if (stddatadir != "") datadir = stddatadir;
    }

    builder->add_from_file(datadir + "/gui.glade"); // May throw exception

    //thr_run(app, builder);
    gui_thread_.reset(new std::thread(&GUI::thr_run, this, app, builder));
}

void GUI::join() {
    gui_thread_->join();
}

void GUI::thr_run(decltype(Gtk::Application::create()) app, decltype(Gtk::Builder::create()) builder) {
    auto window = builder_widget<Gtk::Window>("window1", builder);
    auto vis = builder_widget<Gtk::Viewport>("view_vis", builder);

    // DEBUG:
    addPoint(4, 7, PointType::X);
    addPoint(-4, -7, PointType::X);
    addPoint(0, 3, PointType::X);
    addPoint(9, -2, PointType::X);
    addPoint(7, -2, PointType::X);
    addPoint(5, -2, PointType::X);
    addPoint(5, 1, PointType::X);
    addPoint(-1, 0.09, PointType::X);
    addPoint(2, 2, PointType::X);

    addPoint(0.1, -0.2, PointType::SQUARE);
    addCircle(0.1, -0.2, 2, CircleType::A);
    addCircle(0.1, -0.2, 3, CircleType::B);
    addCircle(0.1, -0.2, 4, CircleType::B);
    addCircle(0.1, -0.2, 5, CircleType::B);
    addCircle(0.1, -0.2, 6, CircleType::B);

    addPoint(-5, -5, PointType::SQUARE);
    addCircle(-5, -5, 1.2, CircleType::A);
    addCircle(-5, -5, 1, CircleType::B);

    GraphArea area(10, 10, -10, -10, *this);

    vis->add(area);
    area.show();

    app->run(*window);
}

unsigned long GUI::addPoint(const double &x, const double &y, const PointType &type) {
    std::lock_guard<std::mutex> lock(mutex_);
    points_.emplace(point_id_next_++, Point_{ .x=x, .y=y, .type=type });
    return point_id_next_-1;
}

bool GUI::removePoint(const unsigned long &id) {
    std::lock_guard<std::mutex> lock(mutex_);
    return points_.erase(id) > 0;
}

void GUI::clearPoints() {
    std::lock_guard<std::mutex> lock(mutex_);
    points_.clear();
}

unsigned long GUI::addCircle(const double &cx, const double &cy, const double &r, const CircleType &type) {
    std::lock_guard<std::mutex> lock(mutex_);
    circles_.emplace(circle_id_next_++, Circle_{ .cx=cx, .cy=cy, .r=r, .type=type });
    return circle_id_next_-1;
}

bool GUI::removeCircle(const unsigned long &id) {
    std::lock_guard<std::mutex> lock(mutex_);
    return circles_.erase(id) > 0;
}

void GUI::clearCircles() {
    std::lock_guard<std::mutex> lock(mutex_);
    circles_.clear();
}

GUI::GraphArea::GraphArea(const double &top, const double &right, const double &bottom, const double &left, GUI &gui) :
    gui_{gui}, bounds_{{top,right,bottom,left}}
{}

bool GUI::GraphArea::on_draw(const Cairo::RefPtr<Cairo::Context> &cr) {

    Gtk::Allocation allocation = get_allocation();
    const int width = allocation.get_width();
    const int height = allocation.get_height();

    // The second of these is typically negative, which is correct because "up" in the graph
    // translates to a lower coordinate on the screen (since positive y coordinates are down the
    // screen).
    const double gwidth = bounds_[RIGHT] - bounds_[LEFT],
          gheight = bounds_[BOTTOM] - bounds_[TOP];

    cr->save();

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
    cr->set_line_width(0.5);
    cr->set_source_rgb(0,0,0);
    cr->move_to(zero_x, 0);
    cr->line_to(zero_x, height);
    cr->stroke();
    cr->move_to(0, zero_y);
    cr->line_to(width, zero_y);
    cr->stroke();

    std::lock_guard<std::mutex> lock(gui_.mutex_);
    for (auto &pt : gui_.points_)
        drawPoint(cr, pt.second, trans);

    for (auto &circ : gui_.circles_)
        drawCircle(cr, circ.second, trans);

    cr->restore();
    return true;
}

void GUI::GraphArea::drawPoint(const Cairo::RefPtr<Cairo::Context> &cr, const Point_ &point, const Cairo::Matrix &trans) {
    // Get the user-space coordinates of the point
    double ux = point.x, uy = point.y;
    trans.transform_point(ux, uy);

    cr->save();
    cr->set_line_width(2.0);

    switch (point.type) {
        case PointType::X:
            {
                // We're drawing 45-degree lines, so get the linear distance on each dimension:
                double edge = gui_.pointMarkerRadius * sqrt(2.0);

                cr->set_source_rgb(1.0, 0, 0); // red

                cr->move_to(ux-0.5*edge, uy-0.5*edge); // Upper left
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
                cr->rectangle(ux-0.5*edge, uy-0.5*edge, edge, edge);
                break;
            }
        case PointType::CROSS:
            {
                cr->set_source_rgb(0.0, 0.25, 1.0); // Blue with a hint of green
                // Horizontal line:
                cr->move_to(ux - gui_.pointMarkerRadius, uy);
                cr->rel_line_to(2*gui_.pointMarkerRadius, 0);

                // Verical line:
                cr->move_to(ux, uy - gui_.pointMarkerRadius);
                cr->rel_line_to(0, 2*gui_.pointMarkerRadius);
                break;
            }
    }

    cr->stroke();
    cr->restore();
}

void GUI::GraphArea::drawCircle(const Cairo::RefPtr<Cairo::Context> &cr, const Circle_ &circle, const Cairo::Matrix &trans) {
    // Apply the full transformation: then we can just draw the circle in *graph* space which,
    // post-transformation, will be exactly the desired oval.
    cr->save();
    cr->transform(trans);

    // Set up the circle in the transformed (graph) coordinates, but *don't* draw it yet (if we do,
    // the line width will also be transformed, which we don't want--it would be too thin on the
    // more compressed axis).
    cr->arc(circle.cx, circle.cy, circle.r, 0.0, 2*M_PI);

    // Undo the transformation, then change the color and line width and actually draw the circle
    cr->restore();
    cr->save();
    cr->set_line_width(2.0);
    switch (circle.type) {
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
