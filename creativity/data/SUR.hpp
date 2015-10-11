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

        /** Returns the number of observations for this SUR object.  Note that this is *not* the
         * number of rows of the X matrix as returned by X(): that will generally be `n() *
         * equations().size()`.
         *
         * Note that if an invalid SUR is constructed with multiple models of different sizes, it is
         * undefined which equation's size will be returned (but such a model will throw an
         * exception when attempting to gather the data to solve the model).
         */
        unsigned int n() const { return eqs_.empty() ? 0 : eqs_.front().depVar()->size(); }

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

        /// Returns the vector of beta values for equation `i`
        Eigen::VectorBlock<const Eigen::VectorXd> beta(unsigned i) const;

        /** Returns the number of variables for equation `i` in this SUR object.  `sur.k(i)` is the
         * number of variables in equation `sur.equations()[i]`.
         */
        const unsigned& k(unsigned i) const;

        /** Returns the degrees of freedom for equation `i`.  `sur.df(i)` is the degrees of
         * freedom value for equation `sur.equations()[i]`.
         */
        int df(unsigned i) const;

        /// Returns the covariance estimate of the beta estimators for equation `i`
        Eigen::Block<const Eigen::MatrixXd> covariance(unsigned i) const;

        /** Returns the standard errors (the square roots of the diagonal of covariance()) of the
         * beta estimates for equation `i`.
         *
         * \throws std::logic_error if the model has not been solved yet by calling `solve()`.
         */
        Eigen::VectorBlock<const Eigen::VectorXd> se(unsigned i) const;

        /** Returns the t-ratios for =0 tests for equation `i`, i.e. `beta(i).array() /
         * se(i).array()`
         *
         * \throws std::logic_error if the model has not been solved yet by calling `solve()`.
         */
        Eigen::VectorBlock<const Eigen::VectorXd> tRatios(unsigned i) const;

        /** Returns the p-values of the t-ratios returned by tRatios(`i`)
         *
         * \throws std::logic_error if the model has not been solved yet by calling `solve()`.
         */
        Eigen::VectorBlock<const Eigen::VectorXd> pValues(unsigned i) const;

        /** Returns \f$s^2\f$, the square of the regression standard error, for equation `i`.
         *
         * \throws std::logic_error if the model has not been solved yet by calling `solve()`.
         */
        const double& s2(unsigned i) const;

        /// Returns the residuals for equation `i`
        Eigen::VectorBlock<const Eigen::VectorXd> residuals(unsigned i) const;

        /// Returns the sum-of-squared residuals for equation `i`
        const double& ssr(unsigned i) const;

        /** Returns the \f$R^2\f$ value for the regression for equation `i`.  This is centered if
         * the model contains a constant, uncentered if it does not.  Note that this constant
         * detection depends on the model's hasConstant() method; a model constructed without a
         * constant but with a SimpleVariable that only takes on a single value will not be
         * considered a constant for the purposes of using centering here.
         */
        const double& Rsq(unsigned i) const;

        /** Returns the y data (without solving the model); `gather()` must have been called either
         * explicitly or by calling `solve()`.  This is the y data for the entire model, that is,
         * the y data for each equation stacked together.
         *
         * \throws std::logic_error if neither `gather()` nor `solve()` has been called.
         */
        const Eigen::VectorXd& y() const;

        /** Returns the sub-vector of y associated with equation `i`.
         *
         * \throws std::logic_error if neither `gather()` nor `solve()` has been called.
         */
        Eigen::VectorBlock<const Eigen::VectorXd> y(unsigned i) const;

        /** Returns the X data (without solving the model); `gather()` must have been called either
         * explicitly or by calling `solve()`.  This is the vector of X data blocks for the entire
         * model, that is, element [i] contains the data for equation i.
         *
         * \throws std::logic_error if neither `gather()` nor `solve()` has been called.
         */
        const std::vector<Eigen::MatrixXd>& X() const;

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
        std::vector<Eigen::MatrixXd> X_;

        /// The model sizes
        std::vector<unsigned> k_;

        /// Model beginning offsets (in beta_, var_beta_, etc.)
        std::vector<unsigned> offset_;

        /// Whether solve() has been called, to populate the below
        bool solved_ = false;

        /// The whole beta vector (all equations)
        Eigen::VectorXd beta_;

        /// The estimated covariance of the beta estimators (all equations)
        Eigen::MatrixXd var_beta_;

        /// The standard errors of the betas (all equations)
        Eigen::VectorXd se_;

        /// t-ratios
        Eigen::VectorXd t_ratios_;

        /// p-values
        Eigen::VectorXd p_values_;

        /// Residuals (all equations)
        Eigen::VectorXd residuals_;

        std::vector<double>
            ssr_, ///< SSR for each equation
            s2_, ///< sigma^2 estimates for each equation
            R2_; ///< R^2 value for each equation

        /// Throws a std::logic_error if the model hasn't been gathered.
        void requireGathered() const { if (!gathered_) throw std::logic_error("Cannot access model data before calling gather()"); }

        /// Throws a std::logic_error if the model hasn't been solved.
        void requireSolved() const { if (!solved_) throw std::logic_error("Cannot obtain model estimates before calling solve()"); }
};


template <class... Eqns>
SUR::SUR(Eqns... eqns) {
    add(std::forward<Eqns>(eqns)...);
}

template <class... MoreEqns>
void SUR::add(const Equation &eq1, MoreEqns... eqns) {
    if (offset_.empty()) offset_.push_back(0);
    else offset_.push_back(offset_.back() + k_.back());
    k_.push_back(eq1.numVars());
    eqs_.push_back(eq1);
    add(std::forward<MoreEqns>(eqns)...);
}

template <class... MoreEqns>
void SUR::add(Equation &&eq1, MoreEqns... eqns) {
    if (offset_.empty()) offset_.push_back(0);
    else offset_.push_back(offset_.back() + k_.back());
    k_.push_back(eq1.numVars());
    eqs_.push_back(std::move(eq1));
    add(std::forward<MoreEqns>(eqns)...);
}

}}
