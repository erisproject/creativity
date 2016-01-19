#include <cairomm/surface.h>

namespace creativity { namespace data { namespace graph {

/** Abstract base class for graph destinations.
 *
 * \sa PDF
 * \sa PNG
 * \sa SVG
 */
class Target {
    public:
        /// Default virtual destructor
        virtual ~Target() = default;

        /// Returns the Cairo::Surface for the current image/page.
        virtual Cairo::RefPtr<Cairo::Surface> surface() = 0;

        /** Called to indicate that the current page/surface is finished and should be saved and a
         * new, blank surface/page created.  Typically this means the surface is written as a new
         * page (PDF) or file (PNG).
         */
        virtual void pageDone() = 0;
};

}}}
