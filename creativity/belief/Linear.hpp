#pragma once
#include <Eigen/Core>
#include <memory>
#include <ostream>

namespace creativity { namespace belief {

/** Base class for a linear model with a natural conjugate, normal-gamma prior. 
 */
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
         * Note that Linear subclasses, so long as they aren't storing additional Eigen types, can
         * rely on their default move constructors.
         */
        Linear(Linear &&move) : Linear(move) {}
        /// Default copy constructor
        Linear(const Linear &copy) = default;
        /// Default copy assignment operator
        Linear& operator=(const Linear &copy) = default;
#endif

        /** Constructs a Linear model of `K` parameters and initializes the various variables (beta,
         * s2, V, n) with highly noninformative priors; specifically this model will initialize
         * parameters with:
         *     beta = 0 vector
         *     s2 = `NONINFORMATIVE_S2` (currently 1e+6)
         *     V = identity matrix
         *     n = `NONINFORMATIVE_N` (currently 1e-3)
         *
         * Unlike calling the full parameter constructor with the above values, this also keeps
         * track of the fact that the model is non-informative, which has two effects:
         *
         * - the first update() call with return a model with `n` set to the number of rows in the
         *   updated data without adding the initial small n value.
         * - noninformative() will return true.
         */
        Linear(unsigned int K);

        // NB: if changing these constants, also change the above constructor documentation
        static constexpr double
            /** The value of `n` for a default noninformative model constructed using
             * `Linear(unsigned int)`.
             */
            //
            NONINFORMATIVE_N = 1e-3,
            /// The value of `s2` for a noninformative model constructed using `Linear(unsigned int)`
            NONINFORMATIVE_S2 = 1e+6;

        /** Constructs a Linear model with the given parameters.  These parameters will be those
         * used for the prior when updating.
         *
         * \param beta the coefficient mean parameters (which, because of restrictions, might not be
         * the actual means).
         * \param s2 the \f$\sigma^2\f$ value of the error term variance.  Typically the \f$\sigma^2\f$ estimate.
         * \param V the model's V matrix (where \f$s^2 V\f$ is the variance matrix of \f$\beta\f$).
         * \param n the number of data points supporting the other values (which can be a
         * non-integer value).
         * \param V_inv A shared pointer to the inverse of `V`, if already calculated.  If the
         * inverse has not already been calculated, it is better to omit this argument: the inverse
         * will be calculated when needed.
         * \param V_chol_L A shared pointer to the `L` matrix of the cholesky decomposition of `V`
         * (where LL' = V).  If the decomposition has not already been calculated, it is better to
         * omit this argument: the decomposition will be calculated when needed.
         *
         * \throws std::runtime_error if any of (`K >= 1`, `V.rows() == V.cols()`, `K == V.rows()`)
         * are not satisfied (where `K` is determined by the number of rows of `beta`).
         */
        Linear(
                const Eigen::Ref<const Eigen::VectorXd> &beta,
                double s2,
                const Eigen::Ref<const Eigen::MatrixXd> &V,
                double n,
                std::shared_ptr<Eigen::MatrixXd> V_inv = nullptr,
                std::shared_ptr<Eigen::MatrixXd> V_chol_L = nullptr
        );

        /** Default constructor: this constructor exists only to allow Linear objects to be default
         * constructed: default constructed objects are models of 0 parameters; such models will
         * throw an std::logic_error exception if any method other than copy or move assignment is
         * attempted on the object.  This primarily exists so that the following is allowed:
         *
         * Linear m;
         * ...
         * m = properly_constructed_model;
         *
         * and similar constructs where the object needs to be default constructed (such as in STL
         * containers).
         */
        Linear() = default;

        /// Virtual destructor
        virtual ~Linear() = default;

        /** Virtual method called during construction to verify the model size.  If this returns a
         * non-zero value, the given parameters (beta, V for the regular constructor, K for the
         * noninformative constructor) must agree with the returned value.  If this returns 0, beta
         * and V must agree with each other, but any model size >= 1 will be accepted.
         */
        virtual unsigned int fixedModelSize() const;

#define NO_EMPTY_MODEL if (K_ == 0) { throw std::logic_error("Cannot use default constructed model object as a model"); }
        /** Accesses the base distribution means value of beta.  Note that this is *not* necessarily
         * the mean of beta and should not be used for prediction; rather it simply returns the
         * distribution parameter value used, which may well not be the mean if any of the beta
         * values have data restrictions.
         */
        const Eigen::VectorXd& beta() const;

        /** Accesses the s2 value of the model. */
        const double& s2() const;

        /** Accesses the n value of the model. */
        const double& n() const;

        /** Accesses the V value of the model. */
        const Eigen::MatrixXd& V() const;

        /** Accesses the inverse of the V value of the model.  If the inverse has not been
         * calculated yet, this calculates and caches the value before returning it.
         */
        const Eigen::MatrixXd& Vinv() const;

        /** Accesses (calculating if not previously calculated) the "L" matrix of the cholesky
         * decomposition of V, where LL' = V.
         */
        const Eigen::MatrixXd& VcholL() const;

        /** Given a row vector of values \f$X^*\f$, predicts \f$y^*\f$ using the current model values.
         * The default implementation provided by this class simply returns the mean \f$X^*
         * \beta\f$ (the mean of the multivariate \f$t\f$ density for an unrestricted,
         * natural conjugate prior); subclasses may override, particularly for models with
         * restrictions or other model parameter distribution assumptions where the parameter means
         * may be quite different than the distributional mean parameters due to the restrictions.
         *
         * \throws std::logic_error if attempting to call predict() on an empty or noninformative
         * model.  (The default implementation would always simply return 0 for noninformative
         * models, which is not a useful prediction).
         */
        virtual double predict(const Eigen::Ref<const Eigen::RowVectorXd> &Xi);

        /** Draws a vector of \f$\beta\f$ values and \f$s^2\f$ values distributed according to the
         * model's structure.  The default implementation simply returns the posterior beta and s^2
         * values, which is suitable for a natural conjugate, normal-gamma distribution without
         * any model restrictions.
         *
         * \returns a const reference to the vector of values.  This same vector is accessible by
         * calling lastDraw().  Note that this vector is reused for subsequent draw() calls and so
         * should be copied into a new vector when storing past draws is necessary.
         *
         * Subclasses overriding this method in should also consider overriding discard().
         */
        virtual const Eigen::VectorXd& draw();

        /** Returns a reference to the vector of \f$\beta\f$ and \f$s^2\f$ values generated by the
         * last call to draw().  If draw() has not yet been called, the vector will be empty.
         */
        const Eigen::VectorXd& lastDraw() const;

        /** This method performs burn-in draws for subclasses that need such a burn-in period.
         * Since the default draw() implementation of this class returns independent draws, the
         * default implementation of this method does nothing; subclasses overriding draw() should
         * also consider overriding this method.
         */
        virtual void discard(unsigned int burn);

        /** The number of parameters of the model, or 0 if this is not a valid model (i.e. a
         * default-constructed model).
         */
        const unsigned int& K() const { return K_; }

        /** Returns true if this model with initialized as a non-informative prior.
         *
         * \sa Linear(unsigned int)
         */
        const bool& noninformative() const;

        /** Overloaded so that a Linear model can be printed nicely with `std::cout << model`.
         */
        friend std::ostream& operator << (std::ostream &os, const Linear &b);

        /** Combines the stored prior parameters with the provided data to create a new Linear
         * object.
         *
         * with a prior.  The new object is
         * suitable as a 
         *
         * \param X the new X data
         * \param y the new y data
         *
         * X and y must have the same number of rows, but it is permitted that the number of rows be
         * less than the number of parameters.  Updating iteratively with a data set broken into
         * small portions will yield the same posterior as updating once with the full data set.
         */
        Linear update(const Eigen::Ref<const Eigen::VectorXd> &y, const Eigen::Ref<const Eigen::MatrixXd> &X) const __attribute__((warn_unused_result));

    protected:
        /** Called during construction to verify that the given parameters are valid.  Subclasses
         * should override to throw an exception if there are problems.  The default implementation
         * in this class does nothing.
         */
        virtual void verifyParameters() const;

        /** This method calls draw() the given number of times, discarding the results.  Unlike
         * discard(), this always calls draw(), while discard() is meant to only call draw() when
         * needed (i.e. when draws are not independent).  This method should generally not be called
         * directly, but instead invoked by subclasses in overridden discard() methods.
         */
        virtual void discardForce(unsigned int burn);

        /// The column vector prior of coefficient means
        Eigen::VectorXd beta_;
        /// The prior value of \f$s^2\f$, the error term variance estimator
        double s2_;
        /** The prior V matrix, which is the \f$(X^\top X)^{-1}\f$ matrix, *not* the \f$s^2(X^\top
         * X)^{-1}\f$ matrix.  This matrix should be symmetric and positive definite.
         */
        Eigen::MatrixXd V_;

        /** The cached inverse of the prior V matrix, which isn't set until/unless needed. */
        mutable std::shared_ptr<Eigen::MatrixXd> V_inv_;

        /** The cached "L" matrix of the cholesky decomposition of V, where LL' = V. */
        mutable std::shared_ptr<Eigen::MatrixXd> V_chol_L_;

        /// The number of data points supporting this model, which need not be an integer.
        double n_;

        /** True if this model was initialized as a non-informative prior.  Any subclasses changing
         * beta/s2/etc. must take care to reset this appropriately.
         *
         * \sa Linear(unsigned int)
         */
        bool noninformative_{false};

        /// The model size.  If 0, this is a default-constructed object which isn't a valid model.
        unsigned int K_{0};

        /** The `K+1` length vector containing the last random draw of \f$\beta\f$ and \f$s^2\f$
         * produced by the draw() method.  Will be an empty vector before the first call to draw().
         */
        Eigen::VectorXd last_draw_;
};

#undef NO_EMPTY_MODEL

}}
