#pragma once
#include <Eigen/Core>
#include <memory>

using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::VectorXd;
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


        // Let Eigen align things if necessary (this will only matter if KK=2 or KK=4)
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    protected:
        /** Constructs a Linear model with the given priors.
         *
         * \param beta_prior the coefficient means
         * \param s2_prior the model \f$s^2\f$
         * \param V_prior the model's V (\f$s^2 V\f$ is the distribution of \f$\beta\f$).
         * \param n_prior the number of data points supporting the other values
         * \param V_prior_inv A shared pointer to the inverse of `V_prior`, if already calculated.
         * This should normally be omitted: the inverse will be calculated when needed.
         *
         * \throws std::runtime_error if any of (`K >= 1`, `V_.rows() == V_.cols()`, `K ==
         * V_.rows()`) are not satisfied (where `K` is determined by the number of rows of
         * `beta_prior`).
         */
        Linear(const VectorKd &beta_prior, const double &s2_prior, const MatrixKd &V_prior, const double &n_prior,
                const std::shared_ptr<MatrixKd> &V_prior_inv = nullptr)
            : K{beta_prior.rows()}, beta_{beta_prior}, s2_{s2_prior}, V_{V_prior}, V_inv_{V_prior_inv}, n_{n_prior}
        {
            // If we're templated as a Dynamic size class, check that the given matrices conform (if
            // not dynamic, this is guaranteed by the fixed constructor parameter sizes).
            if (KK == Eigen::Dynamic) {
                if (K < 1) throw std::runtime_error("Linear model requires at least one parameter");
                if (V_.rows() != V_.cols()) throw std::runtime_error("Linear requires square V_prior matrix");
                if (K != V_.rows()) throw std::runtime_error("Linear requires beta_prior and V_prior of same number of rows");
            }
        }

        /** Uses the stored prior parameters and the provided data to generate a new linear object
         * with posterior parameters
         *
         * \param X the new X data
         * \param y the new y data
         *
         * X and y must have the same number of rows, but it is permitted that the number of rows be
         * less than the number of parameters.  Updating iteratively with a data set broken into
         * small portions will yield the same posterior as updating once with the full data set.
         */
        Linear<KK> update(VectorXd y, MatrixXd X) const {
            if (y.rows() != X.rows())
                throw std::runtime_error("update(y, X) failed: y and X are non-conformable");
            if (X.cols() != K)
                throw std::runtime_error("update(y, X) failed: X has wrong number of columns");

            MatrixXd Xt = X.transpose();
            MatrixXd XtX = Xt * X;
            if (not V_inv_)
                V_inv_ = std::allocate_shared<MatrixKd>(kd_allocator_, V_.colPivHouseholderQr().inverse());

            std::shared_ptr<MatrixKd> V_post_inv = std::allocate_shared<MatrixKd>(kd_allocator_, *V_inv_ + XtX);

            MatrixKd V_post = V_post_inv->colPivHouseholderQr().inverse();

            MatrixKd beta_post = V_post * (*V_inv_ * beta_ + Xt * y);
            double n_post = n_ + X.rows();
            VectorXd residuals = y - X * beta_;
            double s2_post = (n_ * s2_ + residuals.transpose() * (X * V_ * Xt + MatrixXd::Identity(X.rows(), X.rows())) * residuals) / n_post;

            return Linear<KK>(beta_post, s2_post, V_post, n_post, V_post_inv);
        }

        /// The column vector prior of coefficient means
        Matrix<double, Parameters, 1> beta_;
        /// The prior value of \f$s^2\f$, the error term variance estimator
        double s2_;
        /** The prior V matrix, which is the \f$(X^\top X)^{-1}\f$ matrix, *not* the \f$s^2(X^\top
         * X)^{-1}\f$ matrix.  This matrix should be symmetric and positive definite.
         */
        MatrixKd V_;
        /** The inverse of the prior V matrix, which isn't set until/unless needed. */
        std::shared_ptr<MatrixKd> V_inv_;

        /** Allocator for new MatrixKd, used when storing in shared_ptr */
        static Eigen::aligned_allocator<MatrixKd> kd_allocator_;

        /** The prior n, which need not be an integer.
         */
        double n_;
};

}}
