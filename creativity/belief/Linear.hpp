#pragma once
#include <Eigen/Core>
#include <memory>
#include <ostream>

namespace creativity { namespace belief {

/** Base class for a linear model.  This class is simply a container for the common data requires by
 * a linear model plus a simple prediction method that predicts using the current model parameter
 * values.
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
        /** Move assignment operator for Eigen versions before 3.3.  Eigen 3.2 and earlier don't
         * have proper move support, so for Eigen <3.3 we work around it by changing move assignment
         * into copy assignment (which is, of course, much less efficient).
         */
        //Linear& operator=(Linear &&move) { return operator=(move); }
        /// Default copy assignment operator
        Linear& operator=(const Linear &copy) = default;
#endif

        /** Constructs a Linear model with the given parameters.  These parameters will be those
         * used for the prior when updating.
         *
         * \param beta the coefficient mean parameters (which, because of restrictions, might not be
         * the actual means).
         * \param s2 the \f$\sigma^2\f$ value of the error term variance.  Typically the \f$\sigma^2\f$ estimate.
         * \param V the model's V matrix (where \f$s^2 V\f$ is the distribution of \f$\beta\f$).
         * \param n the number of data points supporting the other values (which can be a
         * non-integer value).
         * \param V_inv A shared pointer to the inverse of `V`, if already calculated.  If the
         * inverse has not already been calculated, it is better to omit this argument: the inverse
         * will be calculated when needed.
         *
         * \throws std::runtime_error if any of (`K >= 1`, `V.rows() == V.cols()`, `K == V.rows()`)
         * are not satisfied (where `K` is determined by the number of rows of `beta`).
         */
        Linear(
                const Eigen::Ref<const Eigen::VectorXd> &beta,
                double s2,
                const Eigen::Ref<const Eigen::MatrixXd> &V,
                double n,
                std::shared_ptr<Eigen::MatrixXd> V_inv = nullptr
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

        /// Given a row vector of values, predicts using the current beta_ values.
        double predict(const Eigen::Ref<const Eigen::RowVectorXd> &Xi) const;

        /** The number of parameters of the model, or 0 if this is not a valid model (i.e. a
         * default-constructed model).
         */
        const unsigned int& K() const { return K_; }

        /** Overloaded so that a Linear model can be printed nicely with `std::cout << model`.
         */
        friend std::ostream& operator << (std::ostream &os, const Linear &b);

    protected:
        /** Called during construction to verify that the given parameters are valid.  Subclasses
         * should override to throw an exception if there are problems.  The default implementation
         * in this class does nothing.
         */
        virtual void verifyParameters() const;

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
        Linear update(const Eigen::Ref<const Eigen::VectorXd> &y, const Eigen::Ref<const Eigen::MatrixXd> &X) const;

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

        /// The number of data points supporting this model, which need not be an integer.
        double n_;

    private:
        // The model size.  If 0, this is a default-constructed object which isn't a valid model.
        unsigned int K_{0};
};

#undef NO_EMPTY_MODEL

}}
