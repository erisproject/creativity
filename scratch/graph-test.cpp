#include "creativity/data/graph/Series.hpp"
#include "creativity/data/graph/PDF.hpp"
#include <map>
#include <limits>
#include <iostream>

using namespace creativity::data::graph;

int main() {
    PDF pdf("graph-test.pdf", 6, 4);

#define NaN std::numeric_limits<double>::quiet_NaN()
#define INF std::numeric_limits<double>::infinity()
    std::map<unsigned, double> myline{
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
    std::map<unsigned, std::pair<double, double>> myci{
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
        {18,{-2,28}},
        {19,{17,21}},
        {20,{14,30}},
        {21,{20,22}}
    };

    std::map<unsigned, std::pair<double, double>> mybiggerci(myci);
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

    std::cout << "min=" << min <<", max=" << max << "\n";


    Series s(pdf, myline.begin()->first, myline.rbegin()->first, min, max);
    s.legend_font.set_family("5yearsoldfont");
    s.legend_font.set_size(5 * Pango::SCALE);
    LineStyle linestyle(RGBA(1,0,0));
    s.addLine(myline, linestyle);
    std::string linemarkup = "TEST line OMG <span foreground=\"blue\">OMG</span> <span background=\"red\" font_weight=\"bold\">OMG</span> <b>OMG</b> OMG <span font=\"Serif 3\">OMG <big>OMG <big>OMG <big>OMG <big>OMG <big>OMG <big>OMG <big>OMG <big>OMG <big>OMG</big></big></big></big></big></big></big></big></big></span>";
    s.addLegendItem(linemarkup, FillStyle(Transparent, linestyle), false);

    s.newPage();

    FillStyle cistyle(RGBA(0,0,1,1./3.));
    s.addRegion(myci, cistyle);
    s.addLine(myline, linestyle);
    s.addLegendItem(linemarkup, FillStyle(Transparent, linestyle), false, FillStyle(cistyle.fill_colour, Black));
    s.addLegendItem("90% confidence boundary", cistyle, false);

    s.newPage();
    FillStyle bigcistyle(RGBA(1,0,0,0.25), LineStyle(RGBA(0,0.5,0,0.5)));
    s.addRegion(myci, cistyle);
    s.addRegion(mybiggerci, bigcistyle);
    s.addLine(myline, linestyle);

#define SET_RGBA(X) set_source_rgba(X.red, X.green, X.blue, X.alpha)
    s.addLegendItem(linemarkup, [&](Cairo::RefPtr<Cairo::Context> ctx, const Cairo::Matrix &m) {
            double x = 0, y = 0, w = 1, h = 1;
            m.transform_point(x, y);
            m.transform_distance(w, h);
            ctx->save();
            ctx->move_to(x, y);
            Series::drawRectangle(ctx, w, h, FillStyle(cistyle.fill_colour, Black), true);
            ctx->SET_RGBA(bigcistyle.fill_colour);
            ctx->paint();
            ctx->SET_RGBA(linestyle.colour);
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
            Series::drawRectangle(ctx, w, h, FillStyle(Transparent, Black), true);
            ctx->rectangle(x, y+0.5*h, w, 0.5*h);
            ctx->SET_RGBA(cistyle.fill_colour);
            ctx->fill();
            ctx->SET_RGBA(bigcistyle.fill_colour);
            ctx->paint();
    });

    s.addLegendItem("68% confidence level", bigcistyle, false);

//    s.addLegendItem("transition1", nested_ci, false, std::list<FillStyle>(1, FillStyle(bigcistyle.fill_colour, Black)));

}
