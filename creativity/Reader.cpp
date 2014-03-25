#include "creativity/Reader.hpp"
#include "creativity/Book.hpp"
#include <eris/Random.hpp>
#include <algorithm>

using eris::Position;

namespace creativity {

const std::vector<double> Reader::default_polynomial{{4., -1.}};
const std::vector<double> Reader::default_weights{{1., .75, .5, .25}};

void Reader::weights(std::vector<double> weights) {
    for (size_t i = 0; i < weights.size(); i++) {
        if (weights[i] < 0)
            throw std::domain_error("Invalid weights: weights[" + std::to_string(i) + "] < 0");
        if (i > 0 and weights[i-1] < weights[i])
            throw std::domain_error("Invalid weights: weights[" + std::to_string(i-1) +
                    "] < weights[" + std::to_string(i) + "]");
    }
    weights_ = std::move(weights);
}

void Reader::polynomial(std::vector<double> coef) {
    // 0th coefficient must be positive
    if (coef.size() >= 1 and coef[0] <= 0) throw std::domain_error("Invalid polynomial: coef[0] <= 0");

    // all coefficients must be finite; first and last non-zero coefficients must be negative
    double first = 0.0, last = 0.0;
    for (size_t i = 0; i < coef.size(); i++) {
        if (not std::isfinite(coef[i]))
            throw std::domain_error("Invalid polynomial: coef[" + std::to_string(i) + "] is not finite");

        if (i > 0 and coef[i] != 0) {
            last = coef[i];
            if (first == 0) first = coef[i];
        }
    }

    if (first > 0) throw std::domain_error("Invalid polynomial: first non-zero, non-constant coefficient must be negative.");
    if (last > 0) throw std::domain_error("Invalid polynomial: last non-zero, non-constant coefficient must be negative.");

    double xlast = 0;
    double p = evalPolynomial(0, coef);
    for (const double &x : {1e-6, .001, .01, .1, 1., 10., 100., 1000., 1e6}) {
        double pnext = evalPolynomial(x, coef);
        if (pnext > p)
            throw std::domain_error("Invalid polynomial: polynomial is increasing: f(" + std::to_string(x) + ") > f(" + std::to_string(xlast) + ")");
    }

    polynomial_ = std::move(coef);
}

const std::vector<double>& Reader::polynomial() const {
    return polynomial_;
}

double Reader::evalPolynomial(const double &x) const {
    return evalPolynomial(x, polynomial());
}
double Reader::evalPolynomial(const double &x, const std::vector<double> &polynomial) {
    double p = 0.0;
    double xi = 1.0;
    for (auto &c : polynomial) {
        p += c * xi;
        xi *= x;
    }
    return p;
}

void Reader::interApply() {
    auto &rng = eris::Random::rng();
    // Flip a (weighted) coin to determine creativity:
    if (std::bernoulli_distribution{writer_prob}(rng)) {

        // The book is centered at the reader's position:
        Position bookPos{position()};

        if (writer_book_sd > 0) {
            // Now add some noise to each dimension
            std::normal_distribution<double> book_noise(0, writer_book_sd);
            for (auto &x : bookPos) {
                x += book_noise(rng);
            }
        }

        simulation()->create<Book>(bookPos, sharedSelf());
    }
}

const std::vector<double>& Reader::weights() const {
    return weights_;
}

double Reader::uNthBook(double d, unsigned long n) const {
    if (n >= weights_.size()) return 0;
    return weights_[n] * uBook(d);
}

double Reader::uBook(double d) const {
    if (d < 0) throw std::domain_error("Invalid distance `" + std::to_string(d) +
            "' passed to " + __func__ + "(double): distance cannot be negative.");

    return std::max(0.0, evalPolynomial(d));
}

}
