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
#if !EIGEN_VERSION_AT_LEAST(3,3,0)
        /** Move constructor for Eigen versions before 3.3.  Eigen 3.2 and earlier don't have proper
         * move support, and the implicit ones break things, so we work around this by providing a
         * Move constructor that just calls the implicit copy constructor.  This, of course, means
         * that for old Eigen versions, almost nothing is saved by moving since we actually copy.
         *
         * Eigen 3.3 adds a proper move constructor, and so we don't need this: the default implicit
         * move constructor should work just fine.
         */
        Linear(Linear &&mv) : Linear(mv) {}
        /// Default move constructor
        Linear(const Linear &copy) = default;
#endif

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
