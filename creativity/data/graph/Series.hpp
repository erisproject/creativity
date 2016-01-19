#include <cairomm/surface.h>
#include "creativity/data/graph/style.hpp"
#include <map>

namespace creativity { namespace data { namespace graph {

/** Class to plot graphs of time series values and/or regions. */
class Series {
    public:
        /// No default constructor
        Series() = delete;

        /** Constructs a new series plot that plots to the given Cairo surface.
         *
         * \param surface the Cairo::Surface object where the graph should be written
         * \param tmin the left-hand-side t value
         * \param tmax the right-hand-side t value
         * \param ymin the minimum y value at the bottom of the plot region
         * \param ymax the maximum y value at the top of the plot region
         */
        Series(Cairo::RefPtr<Cairo::Surface> surface, int tmin, double ymin, int tmax, double ymax);

        /** Plots the pairs of coordinates in the given map to the plot.  Keys should be the `t`
         * values.  If any sequential t values are missing or have non-finite values, the line will
         * be broken around the missing values.
         *
         * \param points the map of time to value coordinates to plot
         * \param style the LineStyle; defaults to a 2-pixel wide black line.
         */
        void addLine(const std::map<unsigned, double> &points,
                LineStyle style = LineStyle("black", 2)
                );

        /** Takes a map of time value to double pairs and plots the region between the pair values.
         * If any t values are missing or have one or more non-finite pair values, the region is
         * broken at that t value.
         *
         * \param intervals the map of time to region boundaries.  The order of the pair values does
         * not matter.
         * \param style FillStyle the fill style of the region.  The border of the fill style will
         * be used only for the top and bottom border of the region (i.e. the sides will have no
         * border).  Defaults to a 1/3 opacity blue fill with no (i.e. fully transparent) border.
         */
        void addRegion(const std::map<unsigned, std::pair<double, double>> intervals,
                FillStyle style = FillStyle("rgba(0,0,255,0.3333)"));

    private:
        Cairo::RefPtr<Cairo::Surface> surface_;
        int tmin_, tmax_;
        double ymin_, ymax_;

};

}}}
