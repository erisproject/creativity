#include "creativity/data/graph/PDF.hpp"

namespace creativity { namespace data { namespace graph {

PDF::PDF(std::string filename, double width_inches, double height_inches) :
    // NB: PdfSurface takes points, which are 1/72 of an inch:
    width_{width_inches*72}, height_{height_inches*72}
{
    pdf_ = Cairo::PdfSurface::create(filename, width_, height_);
}

Cairo::RefPtr<Cairo::Surface> PDF::surface() {
    return pdf_;
}

void PDF::newPage() {
    pdf_->show_page();
}

const double& PDF::width() const { return width_; }
const double& PDF::height() const { return height_; }


}}}
