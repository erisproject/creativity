#include "creativity/data/graph/PDF.hpp"

namespace creativity { namespace data { namespace graph {

PDF::PDF(std::string filename, double width_inches, double height_inches) {
    // NB: PdfSurface takes points, which are 1/72 of an inch:
    double scale_x = width_inches*72, scale_y = height_inches*72;
    pdf_ = Cairo::PdfSurface::create(filename, scale_x, scale_y);
    unit_trans_ = Cairo::scaling_matrix(scale_x, scale_y);
}

Cairo::RefPtr<Cairo::Surface> PDF::surface() {
    return pdf_;
}

const Cairo::Matrix& PDF::unitTransformation() const {
    return unit_trans_;
}

void PDF::newPage() {
    pdf_->show_page();
}


}}}
