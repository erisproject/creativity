#pragma once
#include <Eigen/Core>
#include <Eigen/QR>
#include <memory>
#include <ostream>

using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::MatrixBase;
using Eigen::Ref;

namespace creativity { namespace belief {

/** Base class for a linear model.  This class is simply a container for the common data requires by
 * a linear model plus a simple prediction method that predicts using the current model parameter
 * values.
 *
 * \param KK the number of parameters of the model, or Eigen::Dynamic (the default) if the number of
 * parameters is not known but should be inferred from the constructor.
 */
template <int KK = Eigen::Dynamic, typename = typename std::enable_if<KK >= 1 or KK == Eigen::Dynamic>::type>
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
         *
         * Note that Linear<K> subclasses, so long as they aren't storing additional Eigen types,
         * can rely on their default move constructors.
         */
        Linear(Linear &&move) : Linear(move) {}
        /// Default copy constructor
        Linear(const Linear &copy) = default;
        /** Move assignment operator for Eigen versions before 3.3.  Eigen 3.2 and earlier don't
         * have proper move support, so for Eigen <3.3 we work around it by changing move assignment
         * into copy assignment (which is, of course, much less efficient).
         */
        //Linear& operator=(Linear &&move) { return operator=(move); }
        /// Default copy assignment operator
        Linear& operator=(const Linear &copy) = default;
#endif

        /// A column vector type using the KK template value (either fixed or Eigen::Dynamic)
        typedef Matrix<double, KK, 1> VectorKd;
        /// A row vector type using the KK template value (either fixed or Eigen::Dynamic)
        typedef Matrix<double, 1, KK> RowVectorKd;
        /// A square matrix type using the KK template value (either fixed or Eigen::Dynamic)
        typedef Matrix<double, KK, KK> MatrixKd;
        /// A matrix with a dynamic number of rows and KK (either fixed or dynamic) columns
        typedef Matrix<double, Eigen::Dynamic, KK> MatrixXKd;

        /// A typedef to this class (including template parameters)
        typedef Linear<KK> LinearBase;

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
        Linear(
                const Ref<const VectorKd> &beta_prior,
                double s2_prior,
                const Ref<const MatrixKd> &V_prior,
                double n_prior,
                const std::shared_ptr<MatrixKd> &V_prior_inv = nullptr
        ) : beta_{beta_prior}, s2_{s2_prior}, V_{V_prior}, V_inv_{V_prior_inv}, n_{n_prior}, K_{beta_.rows()}
        {
            // If we're templated as a Dynamic size class, check that the given matrices conform (if
            // not dynamic, this is guaranteed by the fixed constructor parameter sizes).
            if (KK == Eigen::Dynamic) {
                if (K() < 1) throw std::runtime_error("Linear model requires at least one parameter");
                if (V_.rows() != V_.cols()) throw std::runtime_error("Linear requires square V_prior matrix");
                if (K() != V_.rows()) throw std::runtime_error("Linear requires beta_prior and V_prior of same number of rows");
                if (V_inv_ and (V_inv_->rows() != V_inv_->cols() or V_inv_->rows() != K()))
                    throw std::runtime_error("Linear constructed with invalid prior inverse");
            }
        }

        /// Default constructor deleted
        Linear() = delete;

        /// Virtual destructor
        virtual ~Linear() = default;

        /** Accesses the value of the beta prior parameter.  Note that this is *not* necessarily the
         * mean and should not be used for prediction; rather it simply returns the distribution
         * parameter value used by the prior.
         */
        const VectorKd& betaPrior() const { return beta_; }

        /// Given a row vector of values, predicts using the current beta_ values.
        double predict(const Eigen::Ref<const RowVectorKd> &Xi) const {
            // FIXME: this is wrong; prediction should be using some sort of sampling
            return Xi * beta_;
        }

        /// The number of parameters of the model
        const long& K() const { return K_; }

        /** Overloaded so that a Linear model can be printed nicely with `std::cout << model`.
         */
        friend std::ostream& operator << (std::ostream &os, const Linear<KK>& b) {
            os << "Linear<" << b.K() << "> model; prior beta:\n" << b.beta_;
            return os;
        }

        // Let Eigen align things if necessary (this will only matter if KK=2 or KK=4)
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    protected:
        /** Called during construction to verify that the given parameters are valid.  Subclasses
         * should override to throw an exception if there are problems.  The default implementation
         * in this class does nothing.
         */
        virtual void verifyParameters() const {}

        /** Uses the stored prior parameters and the provided data to generate a new linear object
         * with posterior parameters.
         *
         * \param X the new X data
         * \param y the new y data
         *
         * X and y must have the same number of rows, but it is permitted that the number of rows be
         * less than the number of parameters.  Updating iteratively with a data set broken into
         * small portions will yield the same posterior as updating once with the full data set.
         */
        Linear<KK> update(const Ref<const VectorXd> &y, const Ref<const MatrixXKd> &X) const {
            if (y.rows() != X.rows())
                throw std::runtime_error("update(y, X) failed: y and X are non-conformable");
            if (X.cols() != K())
                throw std::runtime_error("update(y, X) failed: X has wrong number of columns");

            MatrixXd Xt = X.transpose();
            MatrixXd XtX = Xt * X;
            if (not V_inv_)
                V_inv_ = std::allocate_shared<MatrixKd>(kd_allocator_, V_.colPivHouseholderQr().inverse());

            // Store (V^{-1} + X^\top X) in a shared pointer so that we can pass it along as V^{-1}
            // to the new Linear object so that if *it* needs an inverse, it can just use this
            // instead of having to reinvert to get back to this value.  This saves some
            // computational time but more importantly, gives more accurate numerical results.
            std::shared_ptr<MatrixKd> V_post_inv = std::allocate_shared<MatrixKd>(kd_allocator_,
                    *V_inv_ + XtX
            );

            MatrixKd V_post = V_post_inv->colPivHouseholderQr().inverse();

            VectorKd beta_post = V_post * (*V_inv_ * beta_ + Xt * y);
            double n_post = n_ + X.rows();
            VectorXd residuals = y - X * beta_;
            double s2_post = (n_ * s2_ + residuals.transpose() * (X * V_ * Xt + MatrixXd::Identity(X.rows(), X.rows())) * residuals) / n_post;

            return Linear<KK>(beta_post, s2_post, V_post, n_post, V_post_inv);
        }

        /// The column vector prior of coefficient means
        VectorKd beta_;
        /// The prior value of \f$s^2\f$, the error term variance estimator
        double s2_;
        /** The prior V matrix, which is the \f$(X^\top X)^{-1}\f$ matrix, *not* the \f$s^2(X^\top
         * X)^{-1}\f$ matrix.  This matrix should be symmetric and positive definite.
         */
        MatrixKd V_;

        /** The inverse of the prior V matrix, which isn't set until/unless needed. */
        mutable std::shared_ptr<MatrixKd> V_inv_; // Mutable because we need to create it during const methods,
                                                  // but the creation doesn't actually meaningfully change anything.

        /** The prior n, which need not be an integer.
         */
        double n_;

        /** Eigen allocator for new MatrixKd objects, used when storing in shared_ptr.  This is
         * really only needed for certain fixed values of KK (2 and 4), but doesn't hurt the rest of
         * the time.
         */
        static Eigen::aligned_allocator<MatrixKd> kd_allocator_;

    private:
        long K_;
};

template <int KK, typename Z> Eigen::aligned_allocator<Matrix<double, KK, KK>> Linear<KK, Z>::kd_allocator_{};

}}
