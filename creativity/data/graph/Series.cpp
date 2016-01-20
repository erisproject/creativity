#include "creativity/data/graph/Series.hpp"
#include <cmath>

#include <iostream>

namespace creativity { namespace data { namespace graph {

constexpr double Series::graph_left_, Series::graph_right_, Series::graph_top_, Series::graph_bottom_;

Series::Series(Target &target, int tmin, int tmax, double ymin, double ymax)
    : target_{target}, tmin_{tmin}, tmax_{tmax}, ymin_{ymin}, ymax_{ymax},
    translate_unit_{target.unitTransformation()},
    translate_graph_{
        /* xx */ (graph_right_-graph_left_)/(tmax_-tmin_), /* xy */ 0.0,
        /* xy */ 0.0, /* yy */ (graph_top_-graph_bottom_)/(ymax-ymin),
        /* x0 */ graph_left_ - tmin*(graph_right_-graph_left_)/(tmax_-tmin_),
        /* y0 */ graph_bottom_ - ymin_*(graph_top_-graph_bottom_)/(ymax_-ymin_)}
{

    if (tmin_ >= tmax_ or ymin_ >= ymax_) throw std::invalid_argument("Series error: invalid graph region: graph must have tmin < tmax and ymin < ymax");

    // translate_graph_ right now translates graph points into [0,1]x[0,1] coordinates; we need to also
    // multiply that by the [0,1]x[0,1] -> native surface space for graph translations.
    translate_graph_.multiply(translate_graph_, translate_unit_);
}

void Series::newPage() {
    target_.newPage();
    page_initialized_ = false;
}

void Series::initializePage() {
    if (page_initialized_) return;
    auto ctx = Cairo::Context::create(target_.surface());
    ctx->save();
    ctx->set_source_rgb(1, 1, 1);
    ctx->paint();
    double x = graph_left_, width = graph_right_ - graph_left_, y = graph_top_, height = graph_bottom_-graph_top_;
    translate_unit_.transform_point(x, y);
    translate_unit_.transform_distance(width, height);
    // Fill with the border colour (we'll fill the graph space overtop of it)
    double padding_tb = graph_buffer_tb_ + graph_border_, padding_lr = graph_buffer_lr_ + graph_border_;
    ctx->rectangle(x-padding_lr, y-padding_tb, width+2*padding_lr, height+2*padding_tb);
    ctx->set_source_rgb(0, 0, 0);
    ctx->fill();
    // Fill the graph space over top
    clipToGraph(ctx);
    ctx->set_source_rgb(1, 1, 1);
    ctx->paint();
    ctx->restore();
    page_initialized_ = true;
}

/// Helpful macro for extracting RGBA values in the 4-argument form Cairo wants
#define rgba(X) X.red, X.green, X.blue, X.alpha

void Series::addLine(const std::map<unsigned, double> &points, const LineStyle &style) {
    initializePage();
    auto end = points.upper_bound(tmax_);
    bool have_prev = false, prev_first = false;
    decltype(points.begin()) prev_it;
    auto ctx = Cairo::Context::create(target_.surface());
    ctx->save();
    clipToGraph(ctx);
    ctx->save();
    ctx->set_matrix(translate_graph_);
    for (auto it = points.lower_bound(tmin_); it != end; it++) {
        bool have_curr = std::isfinite(it->second);
        // If the previous point was the first point in its line segment, but doesn't continue to
        // this one (either because this one is excluded, or there are some missing elements), we
        // need to draw its point.
        if (have_prev and prev_first and (not have_curr or it->first != prev_it->first + 1)) {
            ctx->line_to(prev_it->first, prev_it->second);
        }
        prev_first = false;

        // Now look at the current point: if it continues from the previous point, add a line
        // segment between then; if it doesn't continue (either we skipped observations, or the
        // previous point wasn't finite) we start a new line segment at the current point.
        if (have_curr) {
            if (have_prev and it->first == prev_it->first + 1) {
                ctx->line_to(it->first, it->second);
            }
            else {
                // Start a new segment
                ctx->move_to(it->first, it->second);
                prev_first = true;
            }
        }

        have_prev = have_curr;
        prev_it = it;
    }

    // Finally handle the last value (in case it's a point on its own)
    if (have_prev and prev_first)
        ctx->line_to(prev_it->first, prev_it->second);

    ctx->restore(); // Undo transformation
    ctx->set_line_width(style.thickness);
    ctx->set_line_cap(Cairo::LineCap::LINE_CAP_ROUND);
    ctx->set_source_rgba(rgba(style.colour));
    ctx->stroke();
    ctx->restore();
}


void Series::clipToGraph(Cairo::RefPtr<Cairo::Context> ctx) const {
    double x = graph_left_, y = graph_top_, width = graph_right_-graph_left_, height = graph_bottom_-graph_top_;
    translate_unit_.transform_point(x, y);
    translate_unit_.transform_distance(width, height);

    ctx->rectangle(x-graph_buffer_lr_, y-graph_buffer_tb_, width+2*graph_buffer_lr_, height+2*graph_buffer_tb_);
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

void Series::addRegion(const std::map<unsigned, std::pair<double, double>> &intervals, const FillStyle &style) {
    initializePage();
    auto end = intervals.upper_bound(tmax_);

    auto group_start = intervals.end(), prev_it = intervals.end();
    auto ctx = Cairo::Context::create(target_.surface());
    ctx->save();
    clipToGraph(ctx);

    if (style.fill_colour.alpha > 0) {
        ctx->save();
        ctx->set_matrix(translate_graph_);

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
        ctx->set_source_rgba(rgba(style.fill_colour));
        ctx->fill();

        // Now add line segments for any strays
        for (const auto &stray : strays) {
            ctx->move_to(stray.first, stray.second.first);
            ctx->line_to(stray.first, stray.second.second);
        }

        ctx->restore(); // Undo translation (so stroke width with be right)
        ctx->set_line_width(style.stray_thickness);
        ctx->set_source_rgba(rgba(style.fill_colour));
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

