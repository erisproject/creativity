#pragma once
#include <eris/debug.hpp>
#include <Eigen/Core>
#include <memory>
#include <ostream>
#include <vector>

/** This macro is provided to easily provide versions of update() and weaken() that return a proper
 * Derived type (instead of the base class Linear type), so that expressions such as `derived = devired.update(...)`
 * work as expected.  Without this, the above will fail due to there being no conversion from Linear
 * to Derived.
 *
 * This macro will put the methods in "public:" scope, which means "public:" scope will still be in
 * effect after the macro, so either put it in the public: scope, or at the end of some existing
 * scope.
 */
#define CREATIVITY_LINEAR_DERIVED_COMMON_METHODS(Derived) \
    public: \
        /** Updates the model, as done in Linear::update(), but returns an object of this derived type instead of a Linear object. */ \
        [[gnu::warn_unused_result]] Derived update(const Eigen::Ref<const Eigen::VectorXd> &y, const Eigen::Ref<const Eigen::MatrixXd> &X) const & { Derived u(*this); u.updateInPlace(y, X); return u; } \
        /** Updates the model, as done in Linear::update(), but returns an object of this derived type instead of a Linear object. */ \
        [[gnu::warn_unused_result]] Derived update(const Eigen::Ref<const Eigen::VectorXd> &y, const Eigen::Ref<const Eigen::MatrixXd> &X)      && { updateInPlace(y, X); return std::move(*this); } \
        /** Weakens the model, as done in Linear::weaken(), but returns an object of this derived type instead of a Linear object. */ \
        [[gnu::warn_unused_result]] Derived weaken(double s) const & { Derived w(*this); w.weakenInPlace(s); return w; } \
        /** Weakens the model, as done in Linear::weaken(), but returns an object of this derived type instead of a Linear object. */ \
        [[gnu::warn_unused_result]] Derived weaken(double s)      && { weakenInPlace(s); return std::move(*this); }

namespace creativity { namespace belief {

/// Add a MatrixXdR typedef for a row-major MatrixXd
using MatrixXdR = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

/** Base class for a linear model with a natural conjugate, normal-gamma prior.  If deriving from
 * this class, you almost certainly want to use the
 * CREATIVITY_LINEAR_DERIVED_COMMON_METHODS(DerivationName) macro in your class definition.
 */
class Linear {
    public:
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
        /// Default copy constructor
        Linear(const Linear &copy) = default;
        /// Default copy assignment operator
        Linear& operator=(const Linear &copy) = default;

#ifdef EIGEN_HAVE_RVALUE_REFERENCES
        /// Default move constructor
        Linear(Linear &&move) = default;
        /// Default move assignment
        Linear& operator=(Linear &&move) = default;
#else
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
        /** Move assignment for Eigen versions before 3.3: this simply invokes the copy constructor,
         * but is provided so that subclasses still have implicit move constructors.
         */
        Linear& operator=(Linear &&move) { *this = move; return *this; }
#endif

        /** Constructs a Linear model of `K` parameters and initializes the various variables (beta,
         * s2, V, n) with highly noninformative priors; specifically this model will initialize
         * parameters with:
         *     beta = 0 vector
         *     s2 = `NONINFORMATIVE_S2` (currently 1.0)
         *     V = identity matrix times `NONINFORMATIVE_Vc` (currently 1e+8)
         *     n = `NONINFORMATIVE_N` (currently 1e-3)
         *
         * Note, however, that these values are never used: once the model is updated with enough
         * data, the above are determined entirely by the given data without incoporating the above
         * values.
         *
         * If the optional noninf_X and noninf_y arguments are provided and are non-empty matrices,
         * this stores the given data and the model becomes noninformative, but with initial data
         * that will be incorporated the next time data is added.
         *
         * The first update() call with return a model with `n` set to the number of rows in the
         * updated data (plus noninformative data, if any) without adding the initial small n value.
         *
         * noninformative() will return true, and the model cannot be used for prediction until more data
         * is incorporated via update().
         *
         * \param K the number of model parameters
         * \param noninf_X a matrix of `K` columns of data that this model should be loaded with,
         * once enough additional data is added to make \f$X^\top X\f$ invertible.
         * \param noninf_y the y data associated with `noninf_X`
         */
        explicit Linear(
                const unsigned int K,
                const Eigen::Ref<const MatrixXdR> noninf_X = MatrixXdR(),
                const Eigen::Ref<const Eigen::VectorXd> noninf_y = Eigen::VectorXd()
                );

        /** Constructs a Linear model with the given parameters.  These parameters will be those
         * used for the prior when updating.
         *
         * \param beta the coefficient mean parameters (which, because of restrictions, might not be
         * the actual means).
         *
         * \param s2 the \f$\sigma^2\f$ value of the error term variance.  Typically the \f$\sigma^2\f$ estimate.
         *
         * \param V the model's V matrix (where \f$s^2 V\f$ is the variance matrix of \f$\beta\f$).
         * Note: only the lower triangle of the matrix will be used.
         *
         * \param n the number of data points supporting the other values (which can be a
         * non-integer value).
         *
         * \throws std::runtime_error if any of (`K >= 1`, `V.rows() == V.cols()`, `K == V.rows()`)
         * are not satisfied (where `K` is determined by the number of rows of `beta`).
         */
        Linear(
                const Eigen::Ref<const Eigen::VectorXd> beta,
                double s2,
                const Eigen::Ref<const Eigen::MatrixXd> V,
                double n
              );

        /// Virtual destructor
        virtual ~Linear() = default;

        // NB: if changing these constants, also change the single-int, non-informative constructor documentation
        static constexpr double
            /** The value of `n` for a default noninformative model constructed using
             * `Linear(unsigned int)`.
             */
            //
            NONINFORMATIVE_N = 1e-3,
            /// The value of `s2` for a noninformative model constructed using `Linear(unsigned int)`
            NONINFORMATIVE_S2 = 1.0,
            /// The constant for the diagonals of the V matrix for a noninformative model
            NONINFORMATIVE_Vc = 1e+8;


        /** Virtual method called during construction to verify the model size.  If this returns a
         * non-zero value, the given parameters (beta, V for the regular constructor, K for the
         * noninformative constructor) must agree with the returned value.  If this returns 0, beta
         * and V must agree with each other, but any model size >= 1 will be accepted.
         */
        virtual unsigned int fixedModelSize() const;

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

        /** Accesses (calculating if not previous calculated) the inverse of `VcholL()`.  Note that
         * if VcholL() hasn't been calculated yet, this will calculate it. */
        const Eigen::MatrixXd& VcholLinv() const;

        /** Returns the X data that has been added into this model but hasn't yet been used due to
         * it not being sufficiently large and different enough to achieve full column rank.  The
         * data will be scaled appropriately if weaken() has been called, and so does not
         * necessarily reflect the actual data added.
         *
         * \throws std::logic_error if noninformative() would return false.
         */
        const MatrixXdR& noninfXData() const;

        /** Returns the y data associated with noninfXData().  Like that data, this will be scaled
         * if the model has been weakened.
         *
         * \throws std::logic_error if noninformative() would return false.
         */
        const Eigen::VectorXd& noninfYData() const;

        /** Given a row vector of values \f$X^*\f$, predicts \f$y^*\f$ using the current model
         * values.  The default implementation provided by this class simply returns the mean \f$X^*
         * \beta\f$ (the mean of the multivariate \f$t\f$ density for an unrestricted, natural
         * conjugate prior); subclasses may override, particularly for models with restrictions or
         * other model parameter distribution assumptions where the parameter means may have
         * differences from the restricted distribution means due to the restrictions.
         *
         * \throws std::logic_error if attempting to call predict() on an empty or noninformative
         * model.
         */
        virtual double predict(const Eigen::Ref<const Eigen::RowVectorXd> &Xi);

        /** Draws a vector of \f$\beta\f$ values and \f$h^{-1} = \sigma^2\f$ values distributed
         * according to the model's parameters.  The first `K()` values are the drawn \f$\beta\f$
         * values, the last value is the drawn \f$h^{-1}\f$ value.
         *
         * In particular, this uses a gamma distribution to first draw an h value, then uses that h
         * value to draw multivariate normal beta values.  This means the \f$\beta\f$ values will have a
         * multivariate t distribution with mean `beta()`, covariance parameter `s2()*V()`, and
         * degrees of freedom parameter `n()`.
         *
         * \returns a const reference to the vector of values.  This same vector is accessible by
         * calling lastDraw().  Note that this vector is reused for subsequent draw() calls and so
         * should be copied into a new vector when storing past draws is necessary.
         *
         * Subclasses overriding this method in should also consider overriding discard().
         */
        virtual const Eigen::VectorXd& draw();

        /** Draws a multivariate normal with mean \f$\mu\f$ covariance \f$s^2LL^\top\f$ (i.e. takes
         * a constant and a Cholesky decomposition).
         *
         * \param mu the vector means
         * \param L the Cholesky decomposition matrix
         * \param s a standard deviation multiplier for the Cholesky decomposition matrix.  Typically
         * a \f$\sigma\f$ (NOT \f$\sigma^2\f$) value.  If omitted, defaults to 1 (so that you can
         * just pass the Cholesky decomposition of the full covariance matrix).
         *
         * \returns the random multivariate normal vector.
         *
         * \throws std::logic_error if mu and L have non-conforming sizes
         */
        static Eigen::VectorXd multivariateNormal(
                const Eigen::Ref<const Eigen::VectorXd> &mu,
                const Eigen::Ref<const Eigen::MatrixXd> &L,
                double s = 1.0);

        /** Exception class thrown when draw() is unable to produce an admissable draw.  Not thrown
         * by this class (draws never fail) but available for subclass use.
         */
        class draw_failure : public std::runtime_error {
            public:
                /** Constructor.
                 * \param what the exception message.
                 */
                explicit draw_failure(const std::string &what);
                /** Constructor.
                 * \param what the exception message.
                 * \param model the model to append to the exception message.
                 */
                draw_failure(const std::string &what, const Linear &model);
        };


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

        /** Converting a model to a string returns a human-readable model summary.
         */
        virtual operator std::string() const;

        /** Overloaded so that a Linear model can be printed nicely with `std::cout << model`.
         */
        friend std::ostream& operator << (std::ostream &os, const Linear &b);

        /** The display name of the model class to use when converting it to a summary string.
         * Called by the std::string operator when creating a model summary string.  Defaults to
         * "Linear" but subclasses should override.
         */
        virtual std::string display_name() const;

        /** Using the calling object as a prior, uses the provided data to create a new Linear
         * model.
         *
         * If the prior is a noninformative model, and the data (plus any previously incorporated
         * data) will be checked to see if \f$X^\top X\f$ is full-rank: if so, beta, V, s2, and n
         * will be updated using just the data; otherwise, the data will be stored until a
         * subsequent update provides enough (sufficiently varied) data.
         *
         * \param X the new X data
         * \param y the new y data
         *
         * X and y must have the same number of rows, but it is permitted that the number of rows be
         * less than the number of parameters.  Updating iteratively with a data set broken into
         * small portions will yield exactly the same posterior as updating once with the full data
         * set.
         *
         * Calling this method with y and X having no rows at all is permitted: in such a case the
         * returned object is simply a copy of the calling object, but with various state variables
         * (such as last draw) reset.
         */
        [[gnu::warn_unused_result]]
        Linear update(
                const Eigen::Ref<const Eigen::VectorXd> &y,
                const Eigen::Ref<const Eigen::MatrixXd> &X) const &;

        /** Exactly like the above update() method, but optimized for the case where the caller is
         * an rvalue, typically the result of something like:
         *
         *     new = linearmodel.weaken(1.1).update(y, X);
         *
         * The new object is created by moving the current object rather than copying the current
         * object.
         */
        [[gnu::warn_unused_result]]
        Linear update(
                const Eigen::Ref<const Eigen::VectorXd> &y,
                const Eigen::Ref<const Eigen::MatrixXd> &X) &&;

        /** Returns a new Linear object with exactly the same beta, n, and s2 as the caller, but
         * with its `V()` (and related) matrices scaled so that the standard deviation of the
         * parameters is multiplied by given value.  Thus a value of 2 results in a parameter
         * distribution with the same mean, but with double the standard deviation (and 4 times the
         * variance).
         *
         * If the model is noninformative, but has stored data (i.e. the data isn't yet sufficiently
         * large or varied for \f$X^\top X\f$ to be invertible), the stored data itself is
         * "weakened" by scaling both X and y values by the reciprocal of the given scale value.
         * This ensures that the results of the noninformative will be affected in the exact same
         * way as the results had the stored data been sufficient to make the model informative.
         *
         * This is designed to allow using this belief as a prior (via update()), but with a
         * deliberately weakened weight on the prior.
         */
        [[gnu::warn_unused_result]]
        Linear weaken(double stdev_scale) const &;

        /** Just like the above weaken(), but operates on an rvalue reference: the new object that
         * is returned is the result of moving the current object instead of copying the current
         * object into a new object.
         */
        [[gnu::warn_unused_result]]
        Linear weaken(double stdev_scale) &&;

        /** Returns the parameter names.  If not set explicitly, returns a vector of {"0", "1", "2",
         * ...}.
         */
        const std::vector<std::string>& names() const;

        /** Sets the parameter names to the given vector of names.
         *
         * \param names the new names to set.  Must be either empty (to reset names to default), or
         * of length K(), to reset the current names to the given values.
         *
         * \throws std::domain_error if names is of length other than 0 or K().
         */
        void names(const std::vector<std::string> &names);

    protected:

        /** Weakens the current linear model.  This functionality should only be used internally and
         * by subclasses as required for move and copy update methods; weakening should be
         * considered (externally) as a type of construction of a new object.
         */
        virtual void weakenInPlace(double stdev_scale);

        /** Updates the current linear model in place.  This functionality should only be used
         * internally and by subclasses as required for move and copy update methods; updating
         * should be considered (externally) as a type of construction of a new object.
         */
        virtual void updateInPlace(
                const Eigen::Ref<const Eigen::VectorXd> &y,
                const Eigen::Ref<const Eigen::MatrixXd> &X);

        /** Called to reset any internal object state when creating a new derived object.  This is
         * called automatically at the end of weakenInPlace() and updateInPlace() (themselves called
         * when creating a new updated or weakened object), after the new object's various
         * parameters (beta_, n_, etc.) have been updated.
         *
         * The default base implementation clears the last draw value (as returned by draw() and lastDraw()).
         */
        virtual void reset();

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

        mutable std::shared_ptr<Eigen::MatrixXd>
            /// The cached inverse of the prior V matrix, which isn't set until/unless needed.
            V_inv_,

            /// The cached "L" matrix of the cholesky decomposition of V, where LL' = V.
            V_chol_L_,

            /// The cached inverse of the "L" matrix of the cholesky decomposition of V.
            V_chol_L_inv_;

        /// The number of data points supporting this model, which need not be an integer.
        double n_;

        /** Names associated with the beta parameters.  Uses primarily when printing the model.  If
         * unset (or has length != K()), the names are simply "0", "1", "2", ...
         */
        mutable std::shared_ptr<std::vector<std::string>> beta_names_;

        /// True if beta_names_ is simply set to the default (0 ... K-1)
        mutable bool beta_names_default_ = true;

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

    private:
        // Checks that the given matrices conform; called during construction; throws on error.
        void checkLogic();

        /** Do the in-place update using the given X data when we already have an informative model.
         * This is just like updateInPlace, except it is only called after the parameters have been
         * established by having enough data to make \f$X\top X\f$ invertible.  It also doesn't
         * perform various checks on X and y are comforable.
         */
        void updateInPlaceInformative(
                const Eigen::Ref<const Eigen::VectorXd> &y,
                const Eigen::Ref<const Eigen::MatrixXd> &X
                );

        /** The X data that has been loaded into the model but before X has full column rank.  Such
         * a model is still considered noninformative and cannot be used for results until more data
         * is added.
         *
         * If the model is weakened before becoming informative, this X data and the
         * associated y data is scaled by 1/w, where w is the weakening factor, so that \f$(X\^top
         * X)^{-1}\f$ is scaled by \f$w^2\f$, and that y data stays proportional to X.
         */
        std::shared_ptr<MatrixXdR> noninf_X_;
        /// The y data for a non-informative model.  \sa noninf_X_.
        std::shared_ptr<Eigen::VectorXd> noninf_y_;

};

}}
