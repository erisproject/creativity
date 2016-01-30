#include "creativity/data/graph/Series.hpp"
#include "creativity/data/graph/PDF.hpp"
#include <map>
#include <limits>
#include <unordered_map>
#include <iostream>

using namespace creativity::data::graph;

int main() {
    PDF pdf("graph-test.pdf", 6, 4);

#define NaN std::numeric_limits<double>::quiet_NaN()
#define INF std::numeric_limits<double>::infinity()
    std::map<int, double> myline{
        {1,10.},
        {2,20.},
        {3,25.},
        {4,28.},
        {5,25.},
        {6,26.},
        {7,30.},
        {8,0.},
        {9,-10.},
        {11,-11.},
        {12,-5.},
        {14,14},
        {16,16},
        {18,18},
        {19,19},
        {20,20},
        {21,21}
    };
    std::map<int, std::pair<double, double>> myci{
        {1,{5.,15.}},
        {2,{17.,26.}},
        {3,{45.,5.}},
        {4,{22.,29.}},
        {5,{24.,35.}},
        {6,{0,NaN}},
        {7,{15.,32.}},
        {8,{3.,-3.}},
        {9,{-14.,-6}},
        {11,{-10.,-11.5}},
        {12,{12,INF}},
        {14,{13,15}},
        {16,{16,16}},
        {17,{12,14}},
        {18,{-2,28}},
        {19,{17,21}},
        {20,{14,30}},
        {21,{20,22}}
    };

    std::map<int, std::pair<double, double>> mybiggerci(myci);
    for (auto p : myline) {
        auto &ci = mybiggerci.at(p.first);
        ci.first = p.second + (ci.first-p.second)*1.5;
        ci.second = p.second + (ci.second-p.second)*1.5;
    }

    double min = INF, max = -INF;
    for (auto p : myline) {
        if (std::isfinite(p.second)) {
            min = std::min(min, p.second);
            max = std::max(max, p.second);
        }
    }
    for (auto ci : myci) {
        if (std::isfinite(ci.second.first) and std::isfinite(ci.second.second)) {
            min = std::min(min, std::min(ci.second.first, ci.second.second));
            max = std::max(max, std::max(ci.second.first, ci.second.second));
        }
    }
    for (auto ci : mybiggerci) {
        if (std::isfinite(ci.second.first) and std::isfinite(ci.second.second)) {
            min = std::min(min, std::min(ci.second.first, ci.second.second));
            max = std::max(max, std::max(ci.second.first, ci.second.second));
        }
    }


    Series s(pdf, myline.begin()->first, myline.rbegin()->first, min, max, "<b>Hi there!</b>\n<span weight=\"normal\"><small><small>This is <u>my</u> <i>title</i>.</small></small>\n<small><small><small>(If you want your own title, ask nicely.)</small></small></small></span>");
//    s.updateFontFamily("EB Garamond");
    LineStyle linestyle(RGBA(1,0,0));
    s.addLine(myline, linestyle);
    std::string linemarkup = "TEST line OMG <span foreground=\"blue\">OMG</span> <span background=\"red\" font_weight=\"bold\">OMG</span> <b>OMG</b> OMG <span font=\"Serif 3\">OMG <big>OMG <big>OMG <big>OMG <big>OMG <big>OMG <big>OMG <big>OMG <big>OMG <big>OMG</big></big></big></big></big></big></big></big></big></span>";
    s.addLegendItem(linemarkup, FillStyle(Transparent, linestyle), false, Series::default_legend_box_style); //, FillStyle(White, Black));

    s.newPage();
    s.title = "Title 2";
    s.recalcTicks(10, 8, Series::TickEnds::Replace);
    s.graph_style.fill_colour = RGBA(0.95, 0.95, 0.95);
    s.tick_grid_style = LineStyle(White, 1);
    s.legend_style = FillStyle(RGBA(.9,.9,1), Black);

    FillStyle cistyle(RGBA(0,0,1,1./3.));
    s.addRegion(myci, cistyle);
    s.addLine(myline, linestyle);
    s.addLegendItem(linemarkup, FillStyle(Transparent, linestyle), false, FillStyle(cistyle.fill_colour, Black));
    s.addLegendItem("90% confidence boundary", cistyle, false);

    s.legend_graph_space = -1;

    s.newPage();

    FillStyle bigcistyle(RGBA(1,0,0,0.25), LineStyle(RGBA(0,0.5,0,0.5)));
    s.addLegendItem(linemarkup, [&](Cairo::RefPtr<Cairo::Context> ctx, const Cairo::Matrix &m) {
            double x = 0, y = 0, w = 1, h = 1;
            m.transform_point(x, y);
            m.transform_distance(w, h);
            ctx->save();
            ctx->move_to(x, y);
            Series::drawRectangle(ctx, w, h, FillStyle(cistyle.fill_colour, Series::default_legend_box_style.border), true);
            bigcistyle.fill_colour.applyTo(ctx);
            ctx->paint();
            linestyle.colour.applyTo(ctx);
            ctx->set_line_width(linestyle.thickness);
            ctx->move_to(x, y+0.5*h);
            ctx->rel_line_to(x+w, 0);
            ctx->stroke();
    });

    s.addLegendItem("90% confidence boundary", [&](Cairo::RefPtr<Cairo::Context> ctx, const Cairo::Matrix &m) {
            double x = 0, y = 0, w = 1, h = 1;
            m.transform_point(x, y);
            m.transform_distance(w, h);
            ctx->save();
            ctx->move_to(x, y);
            Series::drawRectangle(ctx, w, h, Series::default_legend_box_style, true);
            ctx->rectangle(x, y+0.5*h, w, 0.5*h);
            cistyle.fill_colour.applyTo(ctx);
            ctx->fill();
            bigcistyle.fill_colour.applyTo(ctx);
            ctx->paint();
    });

    s.addLegendItem("68% confidence level", bigcistyle);

    s.recalcTicks(10, 8, Series::TickEnds::Add);
    s.graph_style.fill_colour = White;
    s.tick_grid_style = LineStyle(RGBA(.95,.95,.95), 1);

    s.title = "My title";
    s.x_label = "t";
    bool first = true;
    bool last_was_left = false;
    for (auto lgs : {-1, 0, 5}) {
        s.legend_graph_space = lgs;
        for (auto p : std::list<std::pair<Series::LegendPosition, std::pair<double, double>>>({
                {Series::LegendPosition::Right, {0.,0.}},
                {Series::LegendPosition::Right, {0.,0.25}},
                {Series::LegendPosition::Right, {0.,0.5}},
                {Series::LegendPosition::Right, {0.,0.75}},
                {Series::LegendPosition::Right, {0.,1.}},
                {Series::LegendPosition::Left, {0.,0.}},
                {Series::LegendPosition::Left, {0.,0.25}},
                {Series::LegendPosition::Left, {0.,0.5}},
                {Series::LegendPosition::Left, {0.,0.75}},
                {Series::LegendPosition::Left, {0.,1.}},
                {Series::LegendPosition::Top, {0.,0.}},
                {Series::LegendPosition::Top, {0.25,0.}},
                {Series::LegendPosition::Top, {0.5,0.}},
                {Series::LegendPosition::Top, {0.75,0.}},
                {Series::LegendPosition::Top, {1.,0.}},
                {Series::LegendPosition::Bottom, {0.,0.}},
                {Series::LegendPosition::Bottom, {0.25,0.}},
                {Series::LegendPosition::Bottom, {0.5,0.}},
                {Series::LegendPosition::Bottom, {0.75,0.}},
                {Series::LegendPosition::Bottom, {1.,0.}},
                {Series::LegendPosition::Inside, {1.,0.}}, // TopRight
                {Series::LegendPosition::Inside, {1.,0.5}}, // MidRight
                {Series::LegendPosition::Inside, {1.,1.}}, // BottomRight
                {Series::LegendPosition::Inside, {0.5,0.}}, // TopMid
                {Series::LegendPosition::Inside, {0.5,0.5}}, // Middle
                {Series::LegendPosition::Inside, {0.5,1.}}, // BottomMid
                {Series::LegendPosition::Inside, {0.,0.}}, // TopLeft
                {Series::LegendPosition::Inside, {0.,0.5}}, // MidLeft
                {Series::LegendPosition::Inside, {0.,1.}} // BottomLeft
        })) {
            if (first) first = false;
            else s.newPage();

            last_was_left = s.y_axis_left = !last_was_left;

            if (p.first == Series::LegendPosition::Inside) {
                s.y_label = "net\nutility\n(utility\n- fixed\nincome)";
                s.y_label_rotated = false;
            }
            else {
                s.y_label = "net utility (utility - fixed income)";
                s.y_label_rotated = true;
            }

            s.legend_position = p.first;
            s.legend_rel_x = p.second.first;
            s.legend_rel_y = p.second.second;

            s.addRegion(myci, cistyle);
            s.addRegion(mybiggerci, bigcistyle);
            s.addLine(myline, linestyle);

            s.finishPage();
        }
    }
    //s.legend_top = -100;

//    s.addLegendItem("transition1", nested_ci, false, std::list<FillStyle>(1, FillStyle(bigcistyle.fill_colour, Black)));

}
