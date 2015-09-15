#pragma once
#include "creativity/data/Equation.hpp"
#include <Eigen/Core>
#include <vector>
#include <algorithm>
#include <ostream>
#include <limits>

namespace creativity { namespace data {

/** Class for running a seemingly-unrelated regressions model. */
class SUR {
    public:
        /** Constructs an SUR solver with one or more models.  Given arguments must be accepted by
         * add(...) and are moved (if temporary) or copied.
         *
         * \sa add()
         */
        template <class... Eqns>
        explicit SUR(Eqns... eqns);

        /** Adds one or more equations to the SUR system.  The first equation is copied; subsequent
         * equations are moved (if temporary) or copied via recursive calls to add().
         *
         * If the existing model has been gathered or solved, any gathered or solved values are
         * discarded.
         */
        template <class... MoreEqns>
        void add(const Equation &eq1, MoreEqns... eqns);

        /** Adds one or more equations to the SUR system.  The first equation is moved; subsequent
         * equations are moved (if temporary) or copied via recursive calls to add().
         *
         * If the existing model has been gathered or solved, any gathered or solved values are
         * discarded.
         */
        template <class... MoreEqns>
        void add(Equation &&eq1, MoreEqns... eqns);

        /** Accesses the equations used for this SUR object.
         */
        const std::vector<Equation>& equations() const;

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

        /** Resets any values obtained or calculated by gather() or solve(), forcing the next
         * gather() or solve() call to recalculate.  This method is called implicitly when adding an
         * equation to a previously-gathered or -solved model.
         */
        void clear();

        /// Returns the vector of beta vectors, where element [i] contains the beta values corresponding to equation i.
        const std::vector<Eigen::VectorBlock<Eigen::VectorXd>>& beta();

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

        /** Returns the y data (without solving the model); calls `gather()` if needed.  This is the
         * y data for the entire model, that is, the y data for each equation stacked together.
         */
        const Eigen::VectorXd& y();

        /** Returns the X data (without solving the model); calls `gather()` if needed.  This is the
         * X data for the entire model, that is, the X data for each equation in block diagonal
         * form.
         */
        const Eigen::MatrixXd& X();

        /** Overloaded so that an SUR object can be sent to an output stream; the output consists of
         * the model followed by the model results if the model has been solved; just the model
         * itself if the model hasn't been solved.
         */
        friend std::ostream& operator<<(std::ostream &os, const SUR &ols);

    protected:
        /// Terminating recursive version of add(), not publically callable.
        void add() {}

        /// The model, given during construction
        std::vector<Equation> eqs_;

        /// Whether gather() has been called, to populate y_ and X_
        bool gathered_ = false;

        /// The y vector generated from the model
        Eigen::VectorXd y_;

        /// The X matrix generated from the model
        Eigen::MatrixXd X_;

        /// Whether solve() has been called, to populate the below
        bool solved_ = false;

        /// The whole beta vector
        Eigen::VectorXd beta_full_;

        /// The vector of betas for the individual equations
        std::vector<Eigen::VectorBlock<Eigen::VectorXd>> beta_;

        /// The estimated covariance of the beta estimators
        Eigen::MatrixXd var_beta_;

        /// Residuals
        Eigen::VectorXd residuals_;

        double ssr_ = std::numeric_limits<double>::quiet_NaN(), ///< SSR
               s2_ = std::numeric_limits<double>::quiet_NaN(), ///< sigma^2 estimate
               R2_ = std::numeric_limits<double>::quiet_NaN(); ///< R^2 value

};


template <class... Eqns>
SUR::SUR(Eqns... eqns) {
    add(std::forward<Eqns>(eqns)...);
}

template <class... MoreEqns>
void SUR::add(const Equation &eq1, MoreEqns... eqns) {
    eqs_.push_back(eq1);
    add(std::forward<MoreEqns>(eqns)...);
}

template <class... MoreEqns>
void SUR::add(Equation &&eq1, MoreEqns... eqns) {
    eqs_.push_back(std::move(eq1));
    add(std::forward<MoreEqns>(eqns)...);
}

}}
