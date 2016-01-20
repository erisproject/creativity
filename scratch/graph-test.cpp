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
        {6,NaN},
        {7,30.},
        {8,0.},
        {9,-10.},
        {11,-11.},
        {12,INF},
        {14,14},
        {16,16},
        {18,18},
        {19,19},
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
    s.addLine(myline, RGBA(1,0,0,1));

    s.newPage();

    s.addRegion(myci);
    s.addLine(myline);

    s.newPage();
    s.addRegion(mybiggerci, FillStyle(RGBA(1,0,0,0.25), LineStyle(RGBA(0,0.5,0,0.75))));
    s.addRegion(myci);
    s.addLine(myline);
}
