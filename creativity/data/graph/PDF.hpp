#include "creativity/data/graph/Target.hpp"

namespace creativity { namespace data { namespace graph {

class PDF : public Target {
    public:
        /// Constructs a new PDF target with the given filename and size
        PDF(std::string filename, double width_inches, double height_inches);

        /// Returns the Cairo::PdfSurface for this PDF
        virtual Cairo::RefPtr<Cairo::Surface> surface() override;

        /// Writes the current surface to disk and starts a new page
        virtual void pageDone() override;

    private:
        Cairo::RefPtr<Cairo::PdfSurface> pdf_;

};

}}}
