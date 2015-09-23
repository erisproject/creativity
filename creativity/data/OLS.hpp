#pragma once
#include "creativity/data/Equation.hpp"
#include <Eigen/Core>
#include <ostream>

namespace creativity { namespace data {

/** Class for running a basic OLS regression.
 *
 * Example:
 *
 *     OLS mymodel(Equation(y) + 1 + x1 + x2);
 *     mymodel.solve();
 *     // Do something with mymodel.beta(), etc.
 *     std::cout << mymodel << "\n";
 */
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
        const Equation& model() const { return model_; }

        /** Returns the number of observations for this OLS object.
         */
        unsigned int n() const { return solved_ ? X_.rows() : model_.depVar()->size(); }

        /** Returns the number of variables for this OLS object.
         */
        unsigned int k() const { return solved_ ? X_.cols() : model_.numVars(); }

        /** Returns the degrees of freedom of the model; this is simply `n() - k()`.  Note that this
         * can be negative, but that such a model is certainly unsolvable.
         */
        int df() const { return n() - k(); }

        /** Calculates and stores the final numerical values from the model.  This is called when
         * needed by `solve()`, and does not typically need to be called explicitly.
         *
         * \throws Variable::SizeError if the sizes of the variables in the model do not match.
         */
        void gather();

        /** Attempts to solve the model, if not already done.
         *
         * \throws RankError if X has too few rows, or has non-full-rank columns.
         */
        void solve();

        /** Returns the vector of coefficients (i.e. the beta vector) that solve the model.
         *
         * \throws std::logic_error if the model has not been solved yet by calling `solve()`.
         */
        const Eigen::VectorXd& beta() const;

        /** Returns the covariance estimate of the beta estimators.
         *
         * \throws std::logic_error if the model has not been solved yet by calling `solve()`.
         */
        const Eigen::MatrixXd& covariance() const;

        /** Returns the standard errors (the square roots of the diagonal of covariance()) of the
         * beta() values.
         *
         * \throws std::logic_error if the model has not been solved yet by calling `solve()`.
         */
        const Eigen::VectorXd& se() const;

        /** Returns the t-ratios for =0 tests, i.e. `beta().array() / se().array()`
         *
         * \throws std::logic_error if the model has not been solved yet by calling `solve()`.
         */
        const Eigen::VectorXd& tRatios() const;

        /** Returns the p-values of the t-ratios returned by tRatios()
         *
         * \throws std::logic_error if the model has not been solved yet by calling `solve()`.
         */
        const Eigen::VectorXd& pValues() const;

        /** Returns \f$s^2\f$, the square of the regression standard error.
         *
         * \throws std::logic_error if the model has not been solved yet by calling `solve()`.
         */
        const double& s2() const;

        /** Returns the residuals
         *
         * \throws std::logic_error if the model has not been solved yet by calling `solve()`.
         */
        const Eigen::VectorXd& residuals() const;

        /** Returns the sum-of-squared residuals
         *
         * \throws std::logic_error if the model has not been solved yet by calling `solve()`.
         */
        const double& ssr() const;

        /** Returns the \f$R^2\f$ value for the regression.  This is centered if the model contains
         * a constant, uncentered if it does not.  Note that this constant detection depends on the
         * model's hasConstant() method; a model constructed without a constant but with a
         * SimpleVariable that only takes on a single value will not be considered a constant for
         * the purposes of using centering here.
         *
         * \throws std::logic_error if the model has not been solved yet by calling `solve()`.
         */
        const double& Rsq() const;

        /** The return value of fStat() */
        struct ftest {
            double f; ///< The F value
            unsigned df_numerator; ///< Numerator d.f.
            unsigned df_denominator; ///< Denominator d.f.
            double p; ///< The p-value of the test
        };

        /** Calculates and returns the F-test of all non-constant coefficients in the model being
         * equal to 0.
         *
         * \throws std::logic_error if the model has not been solved yet by calling `solve()`.
         */
        ftest fTest() const;

        /** Returns the y data used to solve the model. `gather()` must have been called either
         * explicitly or by a previous call to `solve()`.
         *
         * \throws std::logic_error if neither `gather()` nor `solve()` has been called.
         */
        const Eigen::VectorXd& y() const;

        /** Returns the X data used to solve the model.  `gather()` must have been called either explicitly or by a previous
         * call to `solve()`.
         *
         * \throws std::logic_error if neither `gather()` nor `solve()` has been called.
         */
        const Eigen::MatrixXd& X() const;

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

        /// Cached standard errors
        Eigen::VectorXd se_;

        /// Cached t ratios
        Eigen::VectorXd t_ratios_;

        /// Cached p-values
        Eigen::VectorXd p_values_;

        /// Residuals
        Eigen::VectorXd residuals_;

        double ssr_ = 0, ///< SSR
               s2_ = 0, ///< sigma^2 estimate
               R2_ = 0; ///< R^2 value

        /// Throws a std::logic_error if the model hasn't been gathered.
        void requireGathered() const { if (!gathered_) throw std::logic_error("Cannot access model data before calling gather()"); }

        /// Throws a std::logic_error if the model hasn't been solved.
        void requireSolved() const { if (!solved_) throw std::logic_error("Cannot obtain model estimates before calling solve()"); }

};

}}
