#pragma once
#include <cairomm/context.h>
#include <eris/noncopyable.hpp>

namespace creativity { namespace data { namespace graph {

/** Abstract base class for graph destinations.
 *
 * \sa PDF
 * \sa PNG
 * \sa SVG
 */
class Target : private eris::noncopyable {
    public:
        /** Creates and returns the Cairo::Surface on which to draw the graph.
         */
        virtual Cairo::RefPtr<Cairo::Surface> surface() = 0;

        /** Returns a Cairo::Matrix that translates [0,1]x[0,1] coordinates into the drawable image
         * surface.  Coordinate (0,0) should map to the top-left corner, (1,1) should map to the
         * bottom-right corner.
         *
         * The default simply returns `Cairo::scaling_matrix(width(), height())`, but subclasses could
         * override.
         */
        virtual Cairo::Matrix unitTransformation() const { return Cairo::scaling_matrix(width(), height()); }

        /// Returns the surface width, in device units
        virtual const double& width() const = 0;

        /// Returns the surface height, in device units
        virtual const double& height() const = 0;

        /** Called to indicate that the current page/surface is finished and should be saved and a
         * new, blank surface/page created.  Typically this means the surface is written as a new
         * page (PDF) or file (PNG).
         */
        virtual void newPage() = 0;
};

}}}
