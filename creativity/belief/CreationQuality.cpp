#include "creativity/belief/CreationQuality.hpp"
#include <eris/Random.hpp>

namespace creativity { namespace belief {

CreationQuality::CreationQuality(double constant, double log_effort, double sigma)
    : beta_{{constant, log_effort}}, sigma_{sigma}
{}

double CreationQuality::draw(const double &effort) {
    std::normal_distribution<double> rnorm(0.0, sigma_);
    return predict(effort) + rnorm(eris::Random::rng());
}

double CreationQuality::predict(const double &effort) {
    if (effort < 0) throw std::domain_error("Invalid negative effort passed to CreationQuality model");
    return beta_[0] + beta_[1] * log(effort);
}

}}
