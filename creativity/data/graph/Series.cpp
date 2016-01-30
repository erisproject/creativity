#include "creativity/data/graph/Series.hpp"
#include <pangomm/init.h>
#include <pangomm/layout.h>
#include <cmath>

namespace creativity { namespace data { namespace graph {

using namespace std::placeholders;

const RectangleStyle Series::default_legend_box_style(White, LineStyle(Black, 0.5));
const LineStyle Series::default_line_style(Black, 1);
const FillStyle Series::default_region_style(RGBA(0,0,1,0.25));


Series::Series(Target &target, int tmin, int tmax, double ymin, double ymax,
        std::string title, std::string x_label, std::string y_label)
    : title{std::move(title)}, x_label{std::move(x_label)},
    y_label{std::move(y_label)}, target_{target}
{
    setExtents(tmin, tmax, ymin, ymax);
    recalcTicks();
}

Series::~Series() {
    finishPage();
}

void Series::setExtents(int tmin, int tmax, double ymin, double ymax, bool retick) {
    if (tmin >= tmax) throw std::invalid_argument("Invalid graph extents: tmin must be less than tmax");
    if (ymin >= ymax) throw std::invalid_argument("Invalid graph extents: ymin must be less than ymax");
    tmin_ = tmin;
    tmax_ = tmax;
    ymin_ = ymin;
    ymax_ = ymax;

    if (retick) recalcTicks();
}

void Series::recalcTicks(unsigned xmax, unsigned ymax, TickEnds end_mode) {
    recalcTicksT(xmax, end_mode);
    recalcTicksY(ymax, end_mode);
}

void Series::recalcTicksT(unsigned max, TickEnds end_mode) {
    if (max < 2) throw std::logic_error("Series::recalcTicksT: max ticks must be >= 2");

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

    while ((unsigned) ((last - first) / t_incr) >= max) {
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
    t_grid.clear();
    t_grid.emplace(first);
    auto grid_last_it = t_grid.emplace(last).first;
    t_ticks.clear();
    auto ticks_end_it = t_ticks.end();
    if (end_mode == TickEnds::Replace) {
        t_ticks.emplace(tmin_);
        ticks_end_it = t_ticks.emplace_hint(ticks_end_it, tmax_);
        first += t_incr;
        last -= t_incr;
    }
    else if (end_mode == TickEnds::Add) {
        if (4*(first - tmin_) >= t_incr) t_ticks.emplace(tmin_);
        if (4*(tmax_ - last)  >= t_incr) ticks_end_it = t_ticks.emplace_hint(ticks_end_it, tmax_);
    }

    for (int i = first; i <= last; i += t_incr) {
        t_ticks.emplace_hint(ticks_end_it, i);
        t_grid.emplace_hint(grid_last_it, i);
    }
}

void Series::recalcTicksY(unsigned max, TickEnds end_mode) {
    // Now do something similar for y, except we can't start from 1 and go up: we start from 1 and
    // then either go down (1 -> 0.5 -> 0.25 -> 0.1 -> ...) *or* up (1 -> 2 -> 5 -> 10 -> ...)
    y_ticks.clear();

    constexpr std::array<double, 5> base{{1, 2, 2.5, 4, 5}};

    unsigned base_i = 0;
    int exp10 = int(std::log10(ymax_ - ymin_));
    // Start out somewhere close (roughly within an order of magnitude):
    double y_incr = base[base_i] * std::pow(10, exp10);
    double last_y_incr = y_incr;

    int first = -max, last = max;
    unsigned last_ticks = 0;
    while (true) {
        if (not std::isfinite(y_incr) or y_incr <= 0)
            throw std::logic_error("recalcTicksY internal error: non-positive increment");
        int try_first = (int) (ymin_ / y_incr);
        if (ymin_ > 0 and std::fmod(ymin_, y_incr) > 0) try_first++;
        int try_last = (int) (ymax_ / y_incr);
        if (ymax_ < 0 and std::fmod(ymax_, y_incr) < 0) try_last--;

        unsigned ticks = try_last - try_first + 1;
        if (last_ticks > max and ticks <= max) {
            // This step moved from too-many-ticks to an okay number, so we're done
            first = try_first;
            last = try_last;
            break;
        }
        else if (last_ticks > 0 and last_ticks <= max and ticks > max) {
            // This step would move from an okay number to too many ticks, so the last step is what
            // we want.
            y_incr = last_y_incr;
            break;
        }
        else {
            last_ticks = ticks;
            first = try_first;
            last = try_last;
            if (last_ticks > max) {
                // Too many ticks: need a bigger tick increment
                if (++base_i == base.size()) {
                    base_i = 0;
                    exp10++;
                }
            }
            else {
                // Too few: need a smaller tick increment
                if (base_i == 0) {
                    base_i = base.size();
                    exp10--;
                }
                base_i--;
            }
            last_y_incr = y_incr;
            y_incr = base[base_i] * std::pow(10, exp10);
        }
    }


    y_grid.clear();
    y_grid.emplace(first*y_incr);
    auto grid_last_it = y_grid.emplace(last*y_incr).first;
    y_ticks.clear();
    auto ticks_end_it = y_ticks.end();
    if (end_mode == TickEnds::Replace) {
        y_ticks.emplace(ymin_);
        ticks_end_it = y_ticks.emplace_hint(ticks_end_it, ymax_);
        first++;
        last--;
    }
    else if (end_mode == TickEnds::Add) {
        if (first*y_incr - ymin_ >= 0.25*y_incr) y_ticks.emplace(ymin_);
        if (ymax_ - last*y_incr >= 0.25*y_incr) ticks_end_it = y_ticks.emplace_hint(ticks_end_it, ymax_);
    }

    for (int i = first; i <= last; i++) {
        y_ticks.emplace_hint(ticks_end_it, i*y_incr);
        y_grid.emplace_hint(grid_last_it, i*y_incr);
    }
}

void Series::updateFonts(const Pango::FontDescription &new_font) {
    for (auto *fptr : {&title_font, &axis_label_font, &tick_font, &legend_font}) {
        fptr->merge(new_font, true);
    }
}

void Series::updateFontFamily(const std::string &new_font_name) {
    Pango::FontDescription font;
    font.set_family(new_font_name);
    updateFonts(font);
}

void Series::newPage() {
    finishPage();
    target_.newPage();
    page_finished_ = false;
}

void Series::finishPage() {
    if (page_finished_) return;

    auto ctx = Cairo::Context::create(target_.surface());
    Cairo::SaveGuard saver(ctx);
    Pango::init();

    // Paint background:
    background_colour.applyTo(ctx);
    ctx->paint();

    const double &total_width = target_.width(), &total_height = target_.height();

    // Distance from the graph border to the surface edge:
    double graph_left = margin_left, graph_right = margin_right, graph_top = margin_top, graph_bottom = margin_bottom;

    double legend_inner_width = 0, legend_inner_height = 0, legend_outer_width = 0, legend_outer_height = 0;
    if (legend_position != LegendPosition::None) {
        // First we need to figure out how much space we need for the legend area
        for (auto &l : legend_items_) {
            auto legend_size = drawLegendItem(ctx, l, true);
            legend_inner_width = std::max(legend_inner_width, legend_size.first);
            if (legend_inner_height > 0) legend_inner_height += legend_spacing;
            legend_inner_height += legend_size.second;
        }
        if (legend_inner_width > 0) {
            legend_outer_width = legend_inner_width + legend_padding_left + legend_padding_right + legend_style.borderLR();
            legend_outer_height = legend_inner_height + legend_padding_top + legend_padding_bottom + legend_style.borderTB();
        }

        // If the legend is outside the graph, it affects one of the edges:
        if (legend_outer_width > 0) {
            if (legend_position == LegendPosition::Right)
                graph_right += legend_graph_space + legend_outer_width;
            if (legend_position == LegendPosition::Left)
                graph_left += legend_graph_space + legend_outer_width;
            if (legend_position == LegendPosition::Top)
                graph_top += legend_graph_space + legend_outer_height;
            if (legend_position == LegendPosition::Bottom)
                graph_bottom += legend_graph_space + legend_outer_height;
        }
    }

    // If there's a title, draw it, and use its size to determine the graph_top edge:
    if (not title.empty()) {
        Black.applyTo(ctx);
        ctx->move_to(0, margin_top);
        auto title_layout = Pango::Layout::create(ctx);
        title_layout->set_font_description(title_font);
        title_layout->set_width(total_width*Pango::SCALE);
        title_layout->set_ellipsize(Pango::EllipsizeMode::ELLIPSIZE_NONE);
        title_layout->set_alignment(Pango::Alignment::ALIGN_CENTER);
        title_layout->set_markup(title);
        int pw, ph;
        title_layout->get_size(pw, ph);
        title_layout->show_in_cairo_context(ctx);
        graph_top += ph/(double)Pango::SCALE + title_padding;
    }

    auto axis_label_layout = Pango::Layout::create(ctx);
    axis_label_layout->set_font_description(axis_label_font);

    auto tick_label_layout = Pango::Layout::create(ctx);
    tick_label_layout->set_font_description(tick_font);

    // If there are horizontal ticks and/or an x-axis label, they determine graph_bottom
    if (not x_label.empty()) {
        axis_label_layout->set_markup(x_label);
        int lw, lh;
        axis_label_layout->get_size(lw, lh);
        // Can't write it yet: it needs to be centered with the graph, and we don't the final graph
        // location yet.
        graph_bottom += lh / (double)Pango::SCALE + axis_label_padding;
    }
    // Tracks how far the last t tick text extends to the right of its tick line center
    double last_t_tick_extends = 0;
    // Track the space immediately below the graph used by the horizontal axis ticks and tick labels
    double htick_space = 0;
    if (not t_ticks.empty()) {
        if (tick_style.length > 0 and tick_style.thickness > 0)
            htick_space += tick_style.length;

        if (tick_font_colour) {
            double max_height = 0;
            int tvw, tvh;
            // Find the largest tick size; they are *probably* all the same, but just in case:
            for (const auto &t : t_ticks) {
                tick_label_layout->set_text(std::to_string(t));
                tick_label_layout->get_size(tvw, tvh);
                max_height = std::max(max_height, tvh / (double)Pango::SCALE);
            }
            if (max_height > 0)
                htick_space += tick_label_space + max_height;

            // Also check the last t tick to make sure it won't go into the right margin: if it
            // does, we'll need to adjust graph_right so that it doesn't.
            last_t_tick_extends = 0.5 * tvw / (double)Pango::SCALE;
        }

        if (htick_space > 0)
            graph_bottom += htick_space;
    }

    // If there are vertical ticks and/or a y-axis label, they determine either graph_left or
    // graph_right
    double &y_axis_graph_side = (y_axis_left ? graph_left : graph_right);
    if (not y_label.empty()) {
        axis_label_layout->set_markup(y_label);
        int lw, lh;
        axis_label_layout->get_size(lw, lh);
        // Can't write it yet: it needs to be centered with the graph, and we don't have the final
        // graph location yet.
        y_axis_graph_side += (y_label_rotated ? lh : lw) / (double)Pango::SCALE + axis_label_padding;
    }
    // Track the space immediately to the right/left of the graph used by the vertical axis ticks and tick labels
    double vtick_space = 0;
    if (not y_ticks.empty()) {
        if (tick_style.length > 0 and tick_style.thickness > 0)
            vtick_space += tick_style.length;

        if (tick_font_colour) {
            double max_width = 0;
            // Find the largest tick size
            for (const auto &y : y_ticks) {
                std::ostringstream oss;
                oss << y;
                tick_label_layout->set_text(oss.str());
                int tvw, tvh;
                tick_label_layout->get_size(tvw, tvh);
                max_width = std::max(max_width, tvw / (double)Pango::SCALE);
            }

            if (max_width > 0)
                vtick_space += tick_label_space + max_width;
        }

        if (vtick_space > 0)
            y_axis_graph_side += vtick_space;
    }

    // Figure out the matrix for converting data points into interior graph coordinates:
    Cairo::Matrix graph_translation(0,0,0,0,0,0);
    bool redo_graph_translation = false;
    do {
        double left = graph_left + graph_style.borderL() + graph_padding_left,
               right = total_width - graph_right - graph_style.borderR() - graph_padding_right,
               top = graph_top + graph_style.borderT() + graph_padding_top,
               bottom = total_height - graph_bottom - graph_style.borderB() - graph_padding_bottom;
        graph_translation.xx = (right - left) / (tmax_ - tmin_);
        graph_translation.yy = (top - bottom) / (ymax_ - ymin_);
        graph_translation.x0 = left - tmin_*(right - left)/(tmax_ - tmin_);
        graph_translation.y0 = bottom - ymin_*(top - bottom)/(ymax_ - ymin_);

        if (redo_graph_translation) {
            // Already redid it; don't need to do it again.
            redo_graph_translation = false;
        }
        else if (last_t_tick_extends > 0) {
            // Now translate the last tick position to see if it extends into the right margin, and if
            // it does, reduce graph_right to compensate.
            double last_tick_x = *t_ticks.rbegin(), dontcare = 0;
            graph_translation.transform_point(last_tick_x, dontcare);
            double extra = total_width - margin_right - last_tick_x - last_t_tick_extends;
            if (extra < 0) {
                graph_right -= extra;
                redo_graph_translation = true;
            }
        }
    } while (redo_graph_translation);

    double graph_outer_width = total_width - graph_right - graph_left, graph_outer_height = total_height - graph_bottom - graph_top;

    // Draw graph area:
    {
        Cairo::SaveGuard saver(ctx);
        ctx->move_to(graph_left, graph_top);
        drawRectangle(ctx, graph_outer_width, graph_outer_height, graph_style, true);
        // We're now clipped inside it
        // If we have tick grid lines, draw them:
        if (tick_grid_style.colour and not (t_grid.empty() and y_grid.empty())) {
            double left = graph_left, right = total_width-graph_right, top = graph_top, bottom = total_height-graph_bottom;
            for (const auto &t : t_grid) {
                double grid_x = t, dontcare = 0;
                graph_translation.transform_point(grid_x, dontcare);
                ctx->move_to(grid_x, top);
                ctx->line_to(grid_x, bottom);
            }
            for (const auto &y : y_grid) {
                double grid_y = y, dontcare = 0;
                graph_translation.transform_point(dontcare, grid_y);
                ctx->move_to(left, grid_y);
                ctx->line_to(right, grid_y);
            }

            tick_grid_style.applyTo(ctx);
            ctx->stroke();
        }
        // Now draw the graph surface (lines, regions, whatever)
        for (const auto &draw : draw_) {
            Cairo::SaveGuard saver(ctx);
            draw(ctx, graph_translation);
        }
        // Done with the drawing components:
        draw_.clear();
    }

    // Horizontal axis ticks:
    if (not t_ticks.empty()) {
        double tick_top = total_height - graph_bottom;
        for (const auto &t : t_ticks) {
            double tick_x = t, dontcare = 0;
            graph_translation.transform_point(tick_x, dontcare);
            if (tick_style.length > 0 and tick_style.thickness > 0) {
                tick_style.applyTo(ctx);
                ctx->move_to(tick_x, tick_top);
                ctx->rel_line_to(0, tick_style.length);
                ctx->stroke();
            }
            if (tick_font_colour) {
                tick_label_layout->set_text(std::to_string(t));
                int tvw, tvh;
                tick_label_layout->get_size(tvw, tvh);
                double width = (double)tvw / Pango::SCALE;
                ctx->move_to(tick_x - 0.5*width, tick_top + tick_style.length + tick_label_space);
                tick_font_colour.applyTo(ctx);
                tick_label_layout->show_in_cairo_context(ctx);
            }
        }
    }

    // Vertical axis ticks:
    if (not y_ticks.empty()) {
        double tick_right = y_axis_left ? graph_left : total_width - graph_right + tick_style.length;
        for (const auto &y : y_ticks) {
            double dontcare = 0, tick_y = y;
            graph_translation.transform_point(dontcare, tick_y);
            if (tick_style.length > 0 and tick_style.thickness > 0) {
                tick_style.applyTo(ctx);
                ctx->move_to(tick_right, tick_y);
                ctx->rel_line_to(-tick_style.length, 0);
                ctx->stroke();
            }
            if (tick_font_colour) {
                std::ostringstream oss;
                oss << y;
                tick_label_layout->set_text(oss.str());
                int tvw, tvh;
                tick_label_layout->get_size(tvw, tvh);
                double width = (double)tvw / Pango::SCALE, height = (double)tvh / Pango::SCALE;
                ctx->move_to(
                        y_axis_left ? tick_right - tick_style.length - tick_label_space - width : tick_right + tick_label_space,
                        tick_y - 0.5*height);
                tick_font_colour.applyTo(ctx);
                tick_label_layout->show_in_cairo_context(ctx);
            }
        }
    }

    // Axis labels:
    Black.applyTo(ctx);
    if (not x_label.empty()) {
        axis_label_layout->set_markup(x_label);
        int lw, lh;
        axis_label_layout->get_size(lw, lh);
        ctx->move_to(
                graph_left + 0.5*(graph_outer_width - lw/(double)Pango::SCALE),
                total_height - graph_bottom + htick_space + axis_label_padding);
        axis_label_layout->show_in_cairo_context(ctx);
    }
    if (not y_label.empty()) {
        Cairo::SaveGuard saver(ctx);
        if (not y_label_rotated) axis_label_layout->set_alignment(y_axis_left ? Pango::Alignment::ALIGN_RIGHT : Pango::Alignment::ALIGN_LEFT);
        axis_label_layout->set_markup(y_label);
        int lw, lh;
        axis_label_layout->get_size(lw, lh);
        if (y_label_rotated) {
            ctx->rotate_degrees(-90);
            ctx->move_to(
                    graph_bottom - total_height + 0.5*(graph_outer_height - lw/(double)Pango::SCALE),
                    y_axis_left
                        ? graph_left - vtick_space - axis_label_padding - lh/(double)Pango::SCALE
                        : total_width - graph_right + vtick_space + axis_label_padding);
        }
        else {
            ctx->move_to(
                    y_axis_left
                        ? graph_left - vtick_space - axis_label_padding - lw/(double)Pango::SCALE
                        : total_width - graph_right + vtick_space + axis_label_padding,
                    graph_top + 0.5*(graph_outer_height - lh/(double)Pango::SCALE));
        }
        axis_label_layout->show_in_cairo_context(ctx);
    }


    // Draw legend area (this has to be down here, after the graph drawing, because it is
    // potentially drawn on top of the graph):
    if (legend_outer_width > 0 and legend_position != LegendPosition::None) {
        Cairo::SaveGuard saver(ctx);
        double legend_x0 = 0, legend_y0 = 0;
        switch (legend_position) {
            case LegendPosition::Inside:
                legend_x0 = (graph_left + graph_style.borderL() + legend_graph_space)
                    + legend_rel_x * (graph_outer_width - 2*legend_graph_space - legend_outer_width - graph_style.borderLR());
                break;
            case LegendPosition::Right:
                legend_x0 = total_width - margin_right - legend_outer_width;
                break;
            case LegendPosition::Left:
                legend_x0 = margin_left;
                break;
            case LegendPosition::Top:
            case LegendPosition::Bottom:
                legend_x0 = graph_left + legend_rel_x * (graph_outer_width - legend_outer_width);
                break;
            case LegendPosition::None:;
        }
        switch (legend_position) {
            case LegendPosition::Inside:
                legend_y0 = (graph_top + graph_style.borderT() + legend_graph_space)
                    + legend_rel_y * (graph_outer_height - 2*legend_graph_space - legend_outer_height - graph_style.borderTB());
                break;
            case LegendPosition::Right:
            case LegendPosition::Left:
                legend_y0 = graph_top + legend_rel_y * (graph_outer_height - legend_outer_height);
                break;
            case LegendPosition::Top:
                legend_y0 = graph_top - legend_graph_space - legend_outer_height;
                break;
            case LegendPosition::Bottom:
                legend_y0 = total_height - margin_bottom - legend_outer_height;
                break;
            case LegendPosition::None:;
        }

        ctx->move_to(legend_x0, legend_y0);

        // Draw the legend border and background, and clip inside the border:
        drawRectangle(ctx, legend_outer_width, legend_outer_height, legend_style, true);

        // Position us inside the padding, where the first legend item goes:
        double legend_x0_inner = legend_x0 + legend_style.border.thickness + legend_padding_left,
               legend_y0_next = legend_y0 + legend_style.border.thickness + legend_padding_top;

        auto it = legend_items_.begin();
        while (it != legend_items_.end()) {
            auto &l = *it;
            ctx->move_to(legend_x0_inner, legend_y0_next);
            auto size = drawLegendItem(ctx, l);
            legend_y0_next += size.second + legend_spacing;
            if (not std::get<0>(l)) // don't preserve this legend item
                it = legend_items_.erase(it);
            else
                it++;
        }
    }

    page_finished_ = true;
}

void Series::drawRectangle(Cairo::RefPtr<Cairo::Context> ctx, double width, double height, const RectangleStyle &style, bool clip) {
    double top, left;
    ctx->get_current_point(left, top);
    if (style.border.thickness > 0 and style.border.colour) {
        style.border.applyTo(ctx);
        if (style.border_top) {
            ctx->move_to(left, top + 0.5*style.border.thickness);
            ctx->rel_line_to(width, 0);
        }
        if (style.border_bottom) {
            ctx->move_to(left, top + height - 0.5*style.border.thickness);
            ctx->rel_line_to(width, 0);
        }
        if (style.border_left) {
            ctx->move_to(left + 0.5*style.border.thickness, top);
            ctx->rel_line_to(0, height);
        }
        if (style.border_right) {
            ctx->move_to(left + width - 0.5*style.border.thickness, top);
            ctx->rel_line_to(0, height);
        }
        ctx->stroke();
    }
    ctx->rectangle(left + style.borderL(), top + style.borderT(),
            width - style.borderLR(), height - style.borderTB());

    if (clip) ctx->clip_preserve();

    style.fill_colour.applyTo(ctx);
    ctx->fill();
}

void Series::addLine(const std::map<int, double> &points, const LineStyle &style) {
    if (page_finished_) throw std::logic_error("Cannot call addLine() on a finished page");
    draw_.emplace_back(std::bind(&Series::drawLine, this, _1, _2, points, style));
}

void Series::drawLine(Cairo::RefPtr<Cairo::Context> ctx, const Cairo::Matrix &translate_graph, const std::map<int, double> &points, const LineStyle &style) const {
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

    style.applyTo(ctx);
    ctx->set_line_cap(Cairo::LineCap::LINE_CAP_ROUND);

    double clip_top, clip_height;
    {
        double x1, x2, clip_bottom;
        ctx->get_clip_extents(x1, clip_top, x2, clip_bottom);
        clip_height = clip_bottom - clip_top;
    }

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
            {
                Cairo::SaveGuard clip_save(ctx);
                ctx->rectangle(xf, clip_top, xb-xf, clip_height);
                ctx->clip();
                {
                    Cairo::SaveGuard trans_save(ctx);
                    ctx->set_matrix(translate_graph);
                    ctx->move_to(line.front().first, line.front().second);
                    bool first = true;
                    for (const auto &l : line) {
                        if (first) first = false;
                        else ctx->line_to(l.first, l.second);
                    }
                }
                ctx->stroke();
            }
        }
    }
}


void Series::addLegendItem(std::string markup, Series::LegendPainterCallback_t &&painter, bool preserve) {
    if (page_finished_) throw std::logic_error("Cannot call addLegendItem() on a finished page");
    legend_items_.emplace_back(std::move(preserve), std::move(markup), std::move(painter));
}
void Series::addLegendItem(std::string markup, const Series::LegendPainterCallback_t &painter, bool preserve) {
    if (page_finished_) throw std::logic_error("Cannot call addLegendItem() on a finished page");
    legend_items_.emplace_back(preserve, markup, painter);
}

std::pair<double, double> Series::drawLegendItem(Cairo::RefPtr<Cairo::Context> ctx, const LegendItem_t &legend, bool dry_run) const {
    Cairo::SaveGuard saver(ctx);
    double x0, y0;
    ctx->get_current_point(x0, y0);
    // We have to handle the text before the box y position, because the size of the text determines
    // where the box ends up.
    Pango::init();
    Black.applyTo(ctx);
    auto pango_layout = Pango::Layout::create(ctx);
    pango_layout->set_font_description(legend_font);
    pango_layout->set_width(legend_text_max_width*Pango::SCALE);
    pango_layout->set_height(legend_text_max_height*Pango::SCALE);
    pango_layout->set_ellipsize(Pango::EllipsizeMode::ELLIPSIZE_END);
    pango_layout->set_markup(std::get<1>(legend));
    double text_width, text_height;
    {
        int pw, ph;
        pango_layout->get_size(pw, ph);
        text_width = pw/(double)Pango::SCALE;
        text_height = ph/(double)Pango::SCALE;
    }

    std::pair<double, double> size(legend_box_width + legend_box_text_gap + text_width, std::max(legend_box_height, text_height));
    if (dry_run) return size;

    double box_y = y0, text_y = y0;
    double y_diff = 0.5*(text_height - legend_box_height);
    if (y_diff < 0) text_y -= y_diff; // The text box is smaller than the legend box, so move the text down to center it
    else box_y += y_diff; // The text is bigger than the box, so the box needs to be moved down to center it with the text

    ctx->move_to(x0 + legend_box_width + legend_box_text_gap, text_y);
    pango_layout->show_in_cairo_context(ctx);

    {
        Cairo::SaveGuard saver(ctx);
        ctx->save();
        // Add clipping region for the box:
        ctx->rectangle(x0, box_y, legend_box_width, legend_box_height);
        ctx->clip();

        // This matrix will convert [0,1]x[0,1] into the legend box:
        std::get<2>(legend)(ctx, Cairo::Matrix(legend_box_width, 0, 0, legend_box_height, x0, box_y));
    }

    return size;
}

void Series::addLegendItem(std::string markup, const FillStyle &lower, bool preserve, const RectangleStyle &bg) {
    addLegendItem(markup, std::bind(&Series::drawSimpleLegendItem, this, _1, _2, std::move(lower), std::move(bg)), preserve);
}

void Series::drawSimpleLegendItem(Cairo::RefPtr<Cairo::Context> ctx, const Cairo::Matrix &to_box, const FillStyle &lower, const RectangleStyle &bg) const {
    double box_x = 0, box_y = 0, box_width = 1, box_height = 1;
    to_box.transform_point(box_x, box_y);
    to_box.transform_distance(box_width, box_height);
    Cairo::SaveGuard saver(ctx);
    ctx->move_to(box_x, box_y);
    // Draw the box and clip it:
    drawRectangle(ctx, box_width, box_height, bg, true);

    // If there's a lower background, paint it
    if (lower.fill_colour) {
        // The x/width/height are too big, but no matter: it's clipped to the space we want anyway
        ctx->rectangle(box_x, box_y + 0.5*box_height, box_width, 0.5*box_height);
        lower.fill_colour.applyTo(ctx);
        ctx->fill();
    }
    // Lower border line
    if (lower.border.thickness > 0 and lower.border.colour) {
        ctx->move_to(box_x, box_y + 0.5*box_height);
        ctx->rel_line_to(box_width, 0);
        lower.border.applyTo(ctx);
        ctx->stroke();
    }
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

void Series::addRegion(const std::map<int, std::pair<double, double>> &intervals, const FillStyle &style) {
    if (page_finished_) throw std::logic_error("Cannot call addRegion() on a finished page");
    draw_.emplace_back(std::bind(&Series::drawRegion, this, _1, _2, intervals, style));
}

void Series::drawRegion(Cairo::RefPtr<Cairo::Context> ctx, const Cairo::Matrix &translate_graph, const std::map<int, std::pair<double, double>> &intervals, const FillStyle &style) const {
    auto end = intervals.upper_bound(tmax_);
    auto group_start = intervals.end(), prev_it = intervals.end();

    if (style.fill_colour) {
        {
            Cairo::SaveGuard saver(ctx);
            ctx->set_matrix(translate_graph);

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
        }

        ctx->set_line_width(stray_thickness);
        style.fill_colour.applyTo(ctx);
        ctx->stroke();
    }

    // Now draw the border lines (if any)
    if (style.border.colour and style.border.thickness > 0) {
        std::map<int, double> max_line, min_line;
        for (auto it = intervals.lower_bound(tmin_); it != end; it++) {
            const auto &a = it->second.first, &b = it->second.second;
            if (std::isfinite(a) and std::isfinite(b)) {
                max_line.emplace_hint(max_line.end(), it->first, std::max(a, b));
                min_line.emplace_hint(min_line.end(), it->first, std::min(a, b));
            }
        }

        drawLine(ctx, translate_graph, min_line, style.border);
        drawLine(ctx, translate_graph, max_line, style.border);
    }
}


}}}

