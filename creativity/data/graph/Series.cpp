#include "creativity/data/graph/Series.hpp"
#include <pangomm/init.h>
#include <pangomm/layout.h>
#include <cmath>

namespace creativity { namespace data { namespace graph {

Series::Series(Target &target, std::string title_markup, int tmin, int tmax, double ymin, double ymax)
    : title_markup{std::move(title_markup)}, target_{target}, tmin_{tmin}, tmax_{tmax}, ymin_{ymin}, ymax_{ymax}
{
    if (tmin_ >= tmax_ or ymin_ >= ymax_) throw std::invalid_argument("Series error: invalid graph region: graph must have tmin < tmax and ymin < ymax");

    recalcTicks();
}

void Series::recalcTicks(unsigned xmax, unsigned ymax, TickEnds end_mode) {
    if (xmax < 2 or ymax < 2) throw std::logic_error("Series::recalcTicks error: both arguments must be >= 2");
    // The increment to consider (will possibly increase this)
    int t_incr = 1;
    // The first and last ticks to be added; first will be >= tmin_, last will be <= tmax_, and both
    // will be multiples of t_incr.
    int first = tmin_, last = tmax_;
    // Tracks the tick increment: we multiply by the repeating cycle (2,2,5/4,8/5,5/4) (starting at
    // the second value) to get the pattern 1 2 2.5 4 5 10 20 25 40 50 100 200 250 400 500 1000 ...
    // (except we skip 2.5 because it's not an integer) which seems a reasonable set of "round"
    // increments (in particular, they are all the multiples which give a nice round number at least
    // every 4 or 5 steps).
    unsigned n = 0;

    while ((unsigned) ((last - first) / t_incr) >= xmax) {
        // Error out before we overflow:
        if (t_incr > std::numeric_limits<int>::max() / 2) throw std::overflow_error("recalcTicks failed: t increment too large");
        ++n %= 5;
        if (t_incr == 2) { n++; t_incr = 4; } // Skip over 2.5
        else if (n < 2) t_incr *= 2;
        else if (n == 3) t_incr = 8 * (t_incr / 5);
        else t_incr = 5 * (t_incr / 4);
        first = tmin_ + t_incr - 1 - (tmin_ + t_incr - 1) % t_incr;
        last = tmax_ - (tmax_ % t_incr);
    }
    t_ticks.clear();
    if (end_mode == TickEnds::Replace) {
        first += t_incr;
        last -= t_incr;
    }
    if ((end_mode == TickEnds::Replace or end_mode == TickEnds::Add) and tmin_ != first)
        t_ticks.emplace(tmin_);

    for (int i = first; i <= last; i += t_incr) {
        t_ticks.emplace_hint(t_ticks.end(), i);
    }

    if ((end_mode == TickEnds::Replace or end_mode == TickEnds::Add) and tmax_ != first)
        t_ticks.emplace_hint(t_ticks.end(), tmax_);


    // Now do something similar for y, except we can't start from 1 and go up: we start from 1 and
    // then either go down (1 -> 0.5 -> 0.25 -> 0.1 -> ...) *or* up (1 -> 2 -> 5 -> 10 -> ...)
    y_ticks.clear();
    y_ticks.emplace(ymin_);
    y_ticks.emplace(ymax_);
}

Cairo::Matrix Series::translateGraph() const {
    double left = graph_left + graph_style.border.thickness + graph_padding_left,
           right = target_.width() - graph_right - graph_style.border.thickness - graph_padding_right,
           top = graph_top + graph_style.border.thickness + graph_padding_top,
           bottom = target_.height() - graph_bottom - graph_style.border.thickness - graph_padding_bottom;

    return Cairo::Matrix(
            (right - left) / (tmax_ - tmin_), // xx
            0, 0, /* xy, yx */
            (top - bottom) / (ymax_ - ymin_), // yy
            left - tmin_*(right - left)/(tmax_ - tmin_), // x0
            bottom - ymin_*(top - bottom)/(ymax_ - ymin_)); // y0
}

void Series::newPage() {
    target_.newPage();
    page_initialized_ = false;
    legend_next_ = 0;
}

void Series::autospaceTitle(std::string new_title) {
    title_markup = std::move(new_title);
    autospaceTitle();
}

void Series::autospaceTitle() {
    if (page_initialized_) throw std::logic_error("autospaceTitle() must be called before adding anything that initializes a page");

    if (title_markup.empty()) {
        graph_top = title_padding_top + title_padding_bottom;
    }
    else {
        Pango::init();
        auto pango_layout = Pango::Layout::create(Cairo::Context::create(target_.surface()));
        pango_layout->set_font_description(title_font);
        pango_layout->set_width((target_.width() - title_padding_left - title_padding_right)*Pango::SCALE);
        pango_layout->set_ellipsize(Pango::EllipsizeMode::ELLIPSIZE_NONE);
        pango_layout->set_alignment(Pango::Alignment::ALIGN_CENTER);
        pango_layout->set_markup(title_markup);
        int pw, ph;
        pango_layout->get_size(pw, ph);
        graph_top = title_padding_top + title_padding_bottom + ph/(double)Pango::SCALE;
    }
}

void Series::initializePage() {
    if (page_initialized_) return;
    auto ctx = Cairo::Context::create(target_.surface());

    ctx->save();

    // Background:
    background_colour.applyTo(ctx);
    ctx->paint();

    // Border:
    double left = graph_left, right = target_.width()-graph_right, top = graph_top, bottom = target_.height()-graph_bottom;
    ctx->move_to(left, top);
    drawRectangle(ctx, right-left, bottom-top, graph_style);

    // Tick marks; these are in the same colour as the graph style border
    auto translate_graph = translateGraph();
    Pango::init();
    auto pango_layout = Pango::Layout::create(ctx);
    // Horizontal axis ticks:
    if (not t_ticks.empty()) {
        double tick_top = bottom;
        pango_layout->set_font_description(tick_font);
        for (const auto &t : t_ticks) {
            double tick_x = t, dontcare = 0;
            translate_graph.transform_point(tick_x, dontcare);
            if (tick_style.length > 0 and tick_style.thickness > 0) {
                tick_style.applyTo(ctx);
                ctx->move_to(tick_x, tick_top);
                ctx->rel_line_to(0, tick_style.length);
                ctx->stroke();
            }
            if (tick_font_colour.alpha > 0) {
                pango_layout->set_text(std::to_string(t));
                int tvw, tvh;
                pango_layout->get_size(tvw, tvh);
                double width = (double)tvw / Pango::SCALE;
                ctx->move_to(tick_x - 0.5*width, tick_top + tick_style.length);
                tick_font_colour.applyTo(ctx);
                pango_layout->show_in_cairo_context(ctx);
            }
        }
    }
    // Vertical axis ticks:
    if (not y_ticks.empty()) {
        double tick_right = graph_left;
        pango_layout->set_font_description(tick_font);
        for (double y : y_ticks) {
            double dontcare = 0, tick_y = y;
            translate_graph.transform_point(dontcare, tick_y);
            if (tick_style.length > 0 and tick_style.thickness > 0) {
                tick_style.applyTo(ctx);
                ctx->move_to(tick_right, tick_y);
                ctx->rel_line_to(-tick_style.length, 0);
                ctx->stroke();
            }
            if (tick_font_colour.alpha > 0) {
                std::ostringstream oss;
                oss << y;
                pango_layout->set_text(oss.str());
                int tvw, tvh;
                pango_layout->get_size(tvw, tvh);
                double width = (double)tvw / Pango::SCALE, height = (double)tvh / Pango::SCALE;
                ctx->move_to(tick_right - tick_style.length - width, tick_y - 0.5*height);
                tick_font_colour.applyTo(ctx);
                pango_layout->show_in_cairo_context(ctx);
            }
        }
    }

    // Draw the title
    if (not title_markup.empty()) {
        double title_max_height = std::max(graph_top - title_padding_top - title_padding_bottom, 0.);

        Pango::init();
        auto pango_layout = Pango::Layout::create(ctx);
        pango_layout->set_font_description(title_font);
        pango_layout->set_width((target_.width() - title_padding_left - title_padding_right)*Pango::SCALE);
        pango_layout->set_height(title_max_height*Pango::SCALE);
        pango_layout->set_ellipsize(Pango::EllipsizeMode::ELLIPSIZE_END);
        pango_layout->set_alignment(Pango::Alignment::ALIGN_CENTER);
        pango_layout->set_markup(title_markup);
        ctx->move_to(title_padding_left, title_padding_top);
        Black.applyTo(ctx);
        pango_layout->show_in_cairo_context(ctx);
    }

    ctx->restore();
    page_initialized_ = true;

    // Restore any saved legend items
    for (const auto &l : legend_) addLegendItem(l.first, l.second, false);
}

void Series::drawRectangle(Cairo::RefPtr<Cairo::Context> ctx, double width, double height, const FillStyle &style, bool clip) {
    double top, left;
    ctx->get_current_point(left, top);
    if (style.border.thickness > 0 and style.border.colour.alpha > 0) {
        ctx->rectangle(left + 0.5*style.border.thickness, top + 0.5*style.border.thickness,
                width - style.border.thickness, height - style.border.thickness);
        style.border.applyTo(ctx);
        ctx->stroke();
    }
    ctx->rectangle(left + style.border.thickness, top + style.border.thickness,
            width - 2*style.border.thickness, height - 2*style.border.thickness);

    if (clip) ctx->clip_preserve();

    style.fill_colour.applyTo(ctx);
    ctx->fill();
}

void Series::addLine(const std::map<int, double> &points, const LineStyle &style) {
    initializePage();
    auto end = points.upper_bound(tmax_);
    bool have_prev = false;
    // List of list of line segments; the inner list will be contiguous
    std::list<std::list<std::pair<int, double>>> segments;
    for (auto it = points.lower_bound(tmin_); it != end; it++) {
        bool have_curr = std::isfinite(it->second);

        // Look at the current point: if it continues from the previous point, add a line segment
        // between them; if it doesn't continue (either we skipped observations, or the previous
        // point wasn't finite) we start a new line segment at the current point.
        if (have_curr) {
            if (not have_prev or it->first != segments.back().back().first + 1) {
                segments.emplace_back();
            }
            segments.back().push_back(*it);
        }

        have_prev = have_curr;
    }

    auto ctx = Cairo::Context::create(target_.surface());
    ctx->save();
    clipToGraph(ctx);
    style.applyTo(ctx);
    ctx->set_line_cap(Cairo::LineCap::LINE_CAP_ROUND);

    double clip_top, clip_height;
    {
        double x1, x2, clip_bottom;
        ctx->get_clip_extents(x1, clip_top, x2, clip_bottom);
        clip_height = clip_bottom - clip_top;
    }

    auto translate_graph = translateGraph();
    for (const auto &line : segments) {
        if (line.size() == 1) {
            // Single point: no clipping to be done, just draw a line that begins and ends at the
            // same point:
            double x = line.front().first, y = line.front().second;
            translate_graph.transform_point(x, y);
            ctx->move_to(x, y);
            ctx->rel_line_to(0, 0);
            ctx->stroke();
        }
        else {
            double ignore = 0;
            // Set up an x clipping region from front to back
            double xf = line.front().first, xb = line.back().first;
            translate_graph.transform_point(xf, ignore);
            translate_graph.transform_point(xb, ignore);
            ctx->save();
            ctx->rectangle(xf, clip_top, xb-xf, clip_height);
            ctx->clip();
            ctx->save();
            ctx->set_matrix(translate_graph);
            ctx->move_to(line.front().first, line.front().second);
            bool first = true;
            for (const auto &l : line) {
                if (first) first = false;
                else ctx->line_to(l.first, l.second);
            }
            ctx->restore(); // Undo translation
            ctx->stroke();
            ctx->restore(); // Undo clipping
        }
    }

    ctx->restore();
}


void Series::addLegendItem(std::string markup, Series::LegendPainterCallback_t &&painter, bool preserve) {
    addLegendItem(markup, painter, false);
    if (preserve) legend_.emplace_back(std::move(markup), std::move(painter));
}
void Series::addLegendItem(std::string markup, const Series::LegendPainterCallback_t &painter, bool preserve) {
    initializePage();

    auto ctx = Cairo::Context::create(target_.surface());

    double left = target_.width() - graph_right + legend_left,
           top = graph_top + legend_top;
    top += legend_next_;

    // We have to handle the text before the box y position, because the size of the text determines
    // where the box ends up.
    double text_x = left + legend_box_width + legend_box_text_gap,
           text_right = target_.width() - legend_right;
    double text_width = text_right - text_x;

    Pango::init();
    auto pango_layout = Pango::Layout::create(ctx);
    pango_layout->set_font_description(legend_font);
    pango_layout->set_width(text_width*Pango::SCALE);
    pango_layout->set_height(legend_text_max_height*Pango::SCALE);
    pango_layout->set_ellipsize(Pango::EllipsizeMode::ELLIPSIZE_END);
    pango_layout->set_markup(markup);
    double actual_height;
    {
        int pw, ph;
        pango_layout->get_size(pw, ph);
        actual_height = ph/(double)Pango::SCALE;
    }
    if (actual_height < legend_box_height) {
        // The text box is smaller than the legend box, so move the text down to center it
        ctx->move_to(text_x, top + 0.5*(legend_box_height - actual_height));
        legend_next_ += legend_box_height + legend_spacing;
    }
    else {
        // The text box is bigger, so `top` needs adjustment to move the box down (to center it with
        // the text)
        ctx->move_to(text_x, top);
        top += 0.5*(actual_height - legend_box_height);
        legend_next_ += actual_height + legend_spacing;
    }
    Black.applyTo(ctx);
    pango_layout->show_in_cairo_context(ctx);

    // The text is done, and the top-left corner of the box is (left,top).

    // This matrix will convert [0,1]x[0,1] into the legend box:
    Cairo::Matrix to_box(legend_box_width, 0, 0, legend_box_height, left, top);

    ctx->save();
    // Add clipping region for the box:
    ctx->rectangle(left, top, legend_box_width, legend_box_height);
    ctx->clip();

    painter(ctx, to_box);

    ctx->restore(); // Undo clipping

    if (preserve) legend_.emplace_back(markup, painter);
}

void Series::addLegendItem(std::string markup, FillStyle lower, bool preserve, FillStyle bg) {
    initializePage();

    addLegendItem(markup, [this,lower,bg](
                Cairo::RefPtr<Cairo::Context> ctx,
                const Cairo::Matrix &to_box) -> void {

        double box_x = 0, box_y = 0, box_width = 1, box_height = 1;
        to_box.transform_point(box_x, box_y);
        to_box.transform_distance(box_width, box_height);
        ctx->save();
        ctx->move_to(box_x, box_y);
        // Draw the box and clip it:
        drawRectangle(ctx, box_width, box_height, bg, true);

        // If there's a lower background, paint it
        if (lower.fill_colour.alpha > 0) {
            // The x/width/height are too big, but no matter: it's clipped to the space we want anyway
            ctx->rectangle(box_x, box_y + 0.5*(box_height + lower.border.thickness),
                    box_width, 0.5*box_height);
            lower.fill_colour.applyTo(ctx);
            ctx->fill();
        }
        // Lower border line
        if (lower.border.thickness > 0 and lower.border.colour.alpha > 0) {
            ctx->move_to(box_x, box_y + 0.5*box_height);
            ctx->rel_line_to(box_width, 0);
            lower.border.colour.applyTo(ctx);
            ctx->stroke();
        }
        ctx->restore(); // Unclip the legend box
    }, preserve);
}

void Series::clipToGraph(Cairo::RefPtr<Cairo::Context> ctx, bool border) const {
    double left = graph_left, top = graph_top, width = target_.width() - graph_right - graph_left, height = target_.height() - graph_bottom - graph_top;
    if (not border) {
        // Don't include the border in the clip region:
        left += graph_style.border.thickness;
        width -= 2*graph_style.border.thickness;
        top += graph_style.border.thickness;
        height -= 2*graph_style.border.thickness;
    }

    ctx->rectangle(left, top, width, height);
    ctx->clip();
}

void Series::boundRegion(
                Cairo::RefPtr<Cairo::Context> ctx,
                std::map<int, std::pair<double, double>>::const_iterator first,
                std::map<int, std::pair<double, double>>::const_iterator last,
                std::map<int, std::pair<double, double>> &strays) {
    if (first == last) {
        strays.insert(*first);
        return;
    }

    auto it = first;
    ctx->move_to(it->first, std::max(it->second.first, it->second.second));
    // First move along the top edge
    for (it++; it != last; it++) {
        ctx->line_to(it->first, std::max(it->second.first, it->second.second));
    }
    ctx->line_to(it->first, std::max(it->second.first, it->second.second));
    // Now move back along the bottom edge
    for (; it != first; it--) {
        ctx->line_to(it->first, std::min(it->second.first, it->second.second));
    }
    ctx->line_to(it->first, std::min(it->second.first, it->second.second));
    ctx->close_path();
}

void Series::addRegion(const std::map<int, std::pair<double, double>> &intervals, const FillStyle &style, double stray_thickness) {
    initializePage();
    auto end = intervals.upper_bound(tmax_);

    auto group_start = intervals.end(), prev_it = intervals.end();
    auto ctx = Cairo::Context::create(target_.surface());
    ctx->save();
    clipToGraph(ctx);

    if (style.fill_colour.alpha > 0) {
        ctx->save();
        ctx->set_matrix(translateGraph());

        std::map<int, std::pair<double, double>> strays;
        for (auto it = intervals.lower_bound(tmin_); it != end; prev_it = it++) {
            bool have_curr = std::isfinite(it->second.first) and std::isfinite(it->second.second);
            // If we don't have a current element or we have a gap, but there is was an active group
            // before, we need to complete the active group region.
            if (group_start != intervals.end() and (not have_curr or it->first != prev_it->first + 1)) {
                boundRegion(ctx, group_start, prev_it, strays);
                group_start = intervals.end();
            }

            // If we don't have an in-progress group but have a good element, start a group
            if (have_curr and group_start == intervals.end()) {
                group_start = it;
            }
        }
        // Enclose the final group:
        if (group_start != intervals.end()) {
            boundRegion(ctx, group_start, prev_it, strays);
        }

        // Fill the groups
        style.fill_colour.applyTo(ctx);
        ctx->fill();

        // Now add line segments for any strays
        for (const auto &stray : strays) {
            ctx->move_to(stray.first, stray.second.first);
            ctx->line_to(stray.first, stray.second.second);
        }

        ctx->restore(); // Undo translation (so stroke width with be right)
        ctx->set_line_width(stray_thickness);
        style.fill_colour.applyTo(ctx);
        ctx->stroke();
    }

    ctx->restore();

    // Now draw the border lines (if any)
    if (style.border.colour.alpha > 0) {
        std::map<int, double> max_line, min_line;
        for (auto it = intervals.lower_bound(tmin_); it != end; it++) {
            const auto &a = it->second.first, &b = it->second.second;
            if (std::isfinite(a) and std::isfinite(b)) {
                max_line.emplace_hint(max_line.end(), it->first, std::max(a, b));
                min_line.emplace_hint(min_line.end(), it->first, std::min(a, b));
            }
        }

        addLine(min_line, style.border);
        addLine(max_line, style.border);
    }
}


}}}

