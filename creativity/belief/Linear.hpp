#pragma once
#include <Eigen/Core>

using Eigen::Matrix;
using Eigen::MatrixBase;

namespace creativity { namespace belief {

/** Base class for a linear model.  This class is simply a container for the common data requires by
 * a linear model plus a simple prediction method that predicts using the current model parameter
 * values.
 *
 * \param Parameters the number of parameters of the model, or Eigen::Dynamic (the default) to
 * if the number of parameters should be inferred from the constructor.
 */
template <int Parameters = Eigen::Dynamic, typename = typename std::enable_if<Parameters >= 1 or Parameters == Eigen::Dynamic>::type>
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
        const long K;
        /// The compile-time number of parameters, or Eigen::Dynamic if not known at compile time
        constexpr static int KK = Parameters;
        /// Default constructor deleted
        Linear() = delete;
        /// Given a row vector of values, predicts using the current beta_ values.
        template<typename Derived>
        double predict(const MatrixBase<Derived> &X) const {
            return X * beta_;
        }

        /// A row vector type using the Parameters template value (either fixed or Eigen::Dynamic)
        typedef Matrix<double, 1, KK> RowVectorKd;
        /// A column vector type using the Parameters template value (either fixed or Eigen::Dynamic)
        typedef Matrix<double, KK, 1> VectorKd;
        /// A square matrix type using the Parameters template value (either fixed or Eigen::Dynamic)
        typedef Matrix<double, KK, KK> MatrixKd;

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    protected:
        /** Stores the given priors.
         *
         * \throws std::runtime_error if any of (`K >= 1`, `V_.rows() == V_.cols()`, `K ==
         * V_.rows()`) are not satisfied (where `K` is determined by the number of rows of
         * `beta_prior`).
         */
        Linear(const VectorKd &beta_prior, const double &s_prior, const MatrixKd &V_prior, const double &n_prior)
            : K{beta_prior.rows()}, beta_{beta_prior}, s_{s_prior}, V_{V_prior}, n_{n_prior}
        {
            // If we're templated as a Dynamic size class, check that the given matrices conform (if
            // not dynamic, this is guaranteed by the fixed constructor parameter sizes).
            if (KK == Eigen::Dynamic) {
                if (K < 1) throw std::runtime_error("Linear model requires at least one parameter");
                if (V_.rows() != V_.cols()) throw std::runtime_error("Linear requires square V_prior matrix");
                if (K != V_.rows()) throw std::runtime_error("Linear requires beta_prior and V_prior of same number of rows");
            }
        }

        /// The column vector prior of coefficients
        Matrix<double, Parameters, 1> beta_;
        /// The prior value of s
        double s_;
        /** The prior V matrix, which is the \f$(X^\top X)^{-1}\f$ matrix, *not* the \f$s^2(X^\top
         * X)^{-1}\f$ matrix.  This matrix should be symmetric and positive definite.
         */
        Matrix<double, Parameters, Parameters> V_;
        /** The prior n, which need not be an integer.
         */
        double n_;
};

}}
