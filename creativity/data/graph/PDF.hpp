#pragma once
#include "creativity/data/graph/Target.hpp"
#include <cairomm/surface.h>

namespace creativity { namespace data { namespace graph {

class PDF : public Target {
    public:
        /// Constructs a new PDF target with the given filename and size
        PDF(std::string filename, double width_inches, double height_inches);

        /// Creates and returns a Cairo::Context that writes to the current page
        virtual Cairo::RefPtr<Cairo::Surface> surface() override;

        /** Returns the transformation matrix that converts [0,1]x[0,1] coordinates into the full
         * surface area.
         */
        virtual const Cairo::Matrix& unitTransformation() const override;

        /// Writes the current surface to disk and starts a new page
        virtual void newPage() override;

    private:
        Cairo::RefPtr<Cairo::PdfSurface> pdf_;
        Cairo::Matrix unit_trans_;

};

}}}
