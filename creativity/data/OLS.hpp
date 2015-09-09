#pragma once
#include "creativity/data/Equation.hpp"
#include <Eigen/Core>
#include <ostream>

namespace creativity { namespace data {

/** Class for running a basic OLS regression. */
class OLS {
    public:
        /// No default constructor
        OLS() = delete;

        /** Constructs an OLS solver with a model.  The model object is copied.
         */
        OLS(const Equation &model);

        /** Constructs an OLS solver, taking over the given model.
         */
        OLS(Equation &&model);

        /** Accesses the model used for this OLS object.
         */
        const Equation& model() const;

        /** Calculates and stores the final numerical values from the model.  This is called when
         * needed, and does not typically need to be called explicitly.
         *
         * \throws Variable::SizeError if the sizes of the variables in the model do not match.
         */
        void gather();

        /** Attempts to solve the model, if not already done.  This method is called automatically
         * by the various other methods when the solution is needed, and rarely needs to be called
         * explicitly.
         *
         * \throws RankError if X has too few rows, or has non-full-rank columns.
         */
        void solve();

        /// Returns the vector of coefficients (i.e. the beta vector) that solve the model.
        const Eigen::VectorXd& beta();

        /// Returns the covariance estimate of the beta estimators
        const Eigen::MatrixXd& covariance();

        /// Returns the regression standard error, \f$^2\f$
        const double& s2();

        /// Returns the residuals
        const Eigen::VectorXd& residuals();

        /// Returns the sum-of-squared residuals
        const double& ssr();

        /** Returns the \f$R^2\f$ value for the regression.  This is centered if the model contains
         * a constant, uncentered if it does not.  Note that this constant detection depends on the
         * model's hasConstant() method; a model constructed without a constant but with a
         * SimpleVariable that only takes on a single value will not be considered a constant for
         * the purposes of using centering here.
         */
        const double& Rsq();

        /// Returns the y data (without solving the model); calls `gather()` if needed.
        const Eigen::VectorXd& y();

        /// Returns the X data (without solving the model); calls `gather()` if needed.
        const Eigen::MatrixXd& X();

        /** Overloaded so that an OLS object can be sent to an output stream; the output consists of
         * the model followed by the model results if the model has been solved; just the model
         * itself if the model hasn't been solved.
         */
        friend std::ostream& operator<<(std::ostream &os, const OLS &ols);
    protected:
        /// The model, given during construction
        Equation model_;

        /// Whether gather() has been called, to populate y_ and X_
        bool gathered_ = false;

        /// The y vector generated from the model
        Eigen::VectorXd y_;

        /// The X matrix generated from the model
        Eigen::MatrixXd X_;

        /// Whether solve() has been called, to populate the below
        bool solved_ = false;

        /// The beta vector
        Eigen::VectorXd beta_;

        /// The estimated covariance of the beta estimators
        Eigen::MatrixXd var_beta_;

        /// Residuals
        Eigen::VectorXd residuals_;

        double ssr_ = 0, ///< SSR
               s2_ = 0, ///< sigma^2 estimate
               R2_ = 0; ///< R^2 value

};

}}
