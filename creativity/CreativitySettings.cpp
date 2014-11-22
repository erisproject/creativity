#include "creativity/CreativitySettings.hpp"
#include <stdexcept>
#include <cmath>

namespace creativity {

double boundaryFromDensity(uint32_t readers, uint32_t dimensions, double density) {
    if (readers == 0) throw std::logic_error("Cannot calculate boundary when readers == 0");
    if (dimensions == 0) throw std::logic_error("Cannot calculate boundary when dimensions == 0");
    if (density <= 0) throw std::logic_error("Cannot calculate boundary when density <= 0");

    const double r_over_d = readers/density;
    // Calculate the boundaries from the density.  Total hypervolume is (2*boundary)^(dimensions),
    // so to achieve `density` we need boundary set as the solution to:
    //     density = readers / ((2*boundary)^(dimensions))
    // which is:
    //     boundary = 1/2 * (readers / density)^(1/dimensions)
    // thus:
    return 0.5 *
        (dimensions == 1 ? r_over_d :
         dimensions == 2 ? std::sqrt(r_over_d) :
         dimensions == 3 ? std::cbrt(r_over_d) :
         std::pow(r_over_d, 1.0/dimensions));
}

double densityFromBoundary(uint32_t readers, uint32_t dimensions, double boundary) {
    if (readers == 0) throw std::logic_error("Cannot calculate density when readers == 0");
    if (dimensions == 0) throw std::logic_error("Cannot calculate density when dimensions == 0");
    if (boundary <= 0) throw std::logic_error("Cannot calculate density when boundary <= 0");

    return readers / std::pow(2*boundary, dimensions);
}

}
