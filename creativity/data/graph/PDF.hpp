#pragma once
#include "creativity/data/graph/Target.hpp"
#include <cairomm/surface.h>

namespace creativity { namespace data { namespace graph {

/** Graph target that writes to pages of a PDF file. */
class PDF : public Target {
    public:
        /// Constructs a new PDF target with the given filename and size
        PDF(std::string filename, double width_inches, double height_inches);

        /// Creates and returns a Cairo::Context that writes to the current page
        virtual Cairo::RefPtr<Cairo::Surface> surface() override;

        /// Returns the surface width
        virtual const double& width() const override;

        /// Returns the surface height
        virtual const double& height() const override;

        /// Writes the current surface to disk and starts a new page
        virtual void newPage() override;

    private:
        Cairo::RefPtr<Cairo::PdfSurface> pdf_;
        double width_, height_;
};

}}}
