#pragma once
#include <Eigen/Core>

using Eigen::Matrix;

namespace creativity { namespace belief {

/** Base class for a linear model.  This class is simply a container for the common data requires by
 * a linear model.
 *
 * \param Parameters the number of parameters of the model.
 */
template <unsigned int Parameters>
class Linear {
    public:
        /// The number of parameters of the model
        constexpr static unsigned int K = Parameters;
        /// Default constructor deleted
        Linear() = delete;
    protected:
        /// Stores the given priors
        Linear(Matrix<double, K, 1> &&beta_prior, double &&s_prior, Matrix<double, K, K> &&V_prior, double &&n_prior)
            : beta_{std::move(beta_prior)}, s_{s_prior}, V_{std::move(V_prior)}, n_{n_prior}
        {}

        /// The column vector prior of coefficients
        Matrix<double, K, 1> beta_;
        /// The prior value of s
        double s_;
        /** The prior V matrix, which is the \f$(X^\top X)^{-1}\f$ matrix, *not* the \f$s^2(X^\top
         * X)^{-1}\f$ matrix.  This matrix should be symmetric and positive definite.
         */
        Matrix<double, K, K> V_;
        /** The prior n, which need not be an integer.
         */
        double n_;
};

}}
