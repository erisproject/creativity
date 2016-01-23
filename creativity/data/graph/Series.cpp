#include "creativity/data/graph/Series.hpp"
#include <pangomm/init.h>
#include <pangomm/layout.h>
#include <cmath>

namespace creativity { namespace data { namespace graph {

Series::Series(Target &target, std::string title_markup, int tmin, int tmax, double ymin, double ymax)
    : title_markup{std::move(title_markup)}, target_{target}, tmin_{tmin}, tmax_{tmax}, ymin_{ymin}, ymax_{ymax}
{
    if (tmin_ >= tmax_ or ymin_ >= ymax_) throw std::invalid_argument("Series error: invalid graph region: graph must have tmin < tmax and ymin < ymax");
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

    // Draw the border:
    double left = graph_left, right = target_.width()-graph_right, top = graph_top, bottom = target_.height()-graph_bottom;
    ctx->move_to(left, top);
    drawRectangle(ctx, right-left, bottom-top, graph_style);
    ctx->restore();
    page_initialized_ = true;

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

void Series::addLine(const std::map<unsigned, double> &points, const LineStyle &style) {
    initializePage();
    auto end = points.upper_bound(tmax_);
    bool have_prev = false;
    // List of list of line segments; the inner list will be contiguous
    std::list<std::list<std::pair<unsigned, double>>> segments;
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
                std::map<unsigned, std::pair<double, double>>::const_iterator first,
                std::map<unsigned, std::pair<double, double>>::const_iterator last,
                std::map<unsigned, std::pair<double, double>> &strays) {
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

        std::map<unsigned, std::pair<double, double>> strays;
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
        std::map<unsigned, double> max_line, min_line;
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

