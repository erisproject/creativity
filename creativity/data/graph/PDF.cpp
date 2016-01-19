#include "creativity/data/graph/PDF.hpp"

namespace creativity { namespace data { namespace graph {

PDF::PDF(std::string filename, double width_inches, double height_inches) {
    // NB: PdfSurface takes points, which are 1/72 of an inch:
    pdf_ = Cairo::PdfSurface::create(filename, width_inches*72, height_inches*72);
}

Cairo::RefPtr<Cairo::Surface> PDF::surface() {
    return pdf_;
}

void PDF::pageDone() {
    pdf_->show_page();
}


}}}
