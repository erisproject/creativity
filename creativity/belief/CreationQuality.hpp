#pragma once
#include "bayes/LinearModel.hpp"

namespace creativity { namespace belief {

class CreationQuality {
    public:
        /** Creates a authorship quality model \f$Q = \beta_1 + \beta_2 \log(effort) + u\f$ with
         * coefficient values as given.
         *
         * \param constant \f$\beta_1\f$, the constant.
         * \param log_effort \f$\beta_2\f$, the coefficient on \f$\log(effort)\f$.
         * \param sigma the standard deviation of the error term.  \f$u\f$ is normally distributed
         * with mean 0 and this standard deviation.
         */
        CreationQuality(double constant, double log_effort, double sigma);

        /** Draws a quality value given a value of effort, including the random error term.
         * Internally, this method calls `predict(effort)` and adds a random error draw.
         *
         * \throws std::domain_error if effort is less than 0.
         */
        double draw(const double &effort);

        /** Predicts a quality value given a value of effort.  This is similar to draw, but includes
         * no error term.
         *
         * \throws std::domain_error if effort is less than 0.
         */
        double predict(const double &effort);

    private:
        std::array<double, 2> beta_;
        double sigma_;
};

}}
