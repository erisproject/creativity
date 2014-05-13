#pragma once
#include <creativity/belief/Linear.hpp>

namespace creativity { namespace belief {

/** This class represents an author's belief about the lifetime profitability of a work.  The model
 * is of the form:
 *
 * \f$\Pi_b = \beta_1 + \beta_2 q_b^D + \beta_3 firstBook + \beta_4 previousBooks + \beta_5 marketBooks + u\f$
 * where:
 * - \f$Pi_b\f$ is the lifetime profits of the book
 * - \f$q_b\f$ is the quality of the book, which must be non-negative
 * - \f$D\f$ is the dimensionality of the modelled world (e.g. 2 for a two-dimensional world).
 *   \f$q_b\f$ is raised to the dimensionality because changes in it affect the radius of potential
 *   customers, with total customers being proportional to the radius raised to \f$D\f$.
 * - \f$firstBook\f$ is a dummy: 1 if this is the creator's first work, 0 if the creator has other
 *   works.
 * - \f$previousBooks\f$ is the number of previous books the author has created.
 * - \f$marketBooks\f$ is the number of books on the market in the previous period.
 *
 * The following restrictions are imposed on beliefs:
 * - \f$\beta_2 \geq 0\f$ (profit increases with quality)
 * - \f$\beta_5 \leq 0\f$ (more competition means lower profit)
 *
 * These constraints are combined with a natural conjugate prior for the purposes of updating the
 * beliefs via Bayesian econometrics.
 */
class Profit : public Linear<5> {
    public:
        /** Constructs a demand model with the given prior information.
         *
         * \param D the dimensionality of the world.
         * \param beta_prior the prior of the mean values of beta.  Must be a 5-value (column)
         * vector.  Values are in the order given in this class's description.
         * \param s2_prior the prior of \f$s^2\f$ (typically an estimate thereof).
         * \param V_prior the prior covariance matrix of the estimators, *without* premultiplication
         * by \f$\sigma^2\f$.  That is, for a prior from OLS, this is the matrix \f$(X^\top
         * X)^{-1}\f$, not the matrix \f$s^2(X^\top X)^{-1}\f$.  This matrix should be symmetric,
         * positive definite.
         * \param n_prior the number of data points supporting the prior (which need not be an
         * integer).
         */
        Profit(
                const unsigned int &D,
                const VectorKd &beta_prior,
                const double &s2_prior,
                const MatrixKd &V_prior,
                const double &n_prior
              );

        /** Given a set of model parameters, this returns an expected value \f$\Pi_b\f$, the
         * lifetime profit of the book.
         *
         * \param q the quality of the book
         * \param previousBooks the number of previous books created by this book's author.  This
         * parameter also determines the `firstBook` dummy (`= 1` iff `previousBooks == 0`).
         * \param marketBooks the number of books on the market last period
         */
        double predict(const double &q, const unsigned long &previousBooks, const unsigned long &marketBooks) const;

        /** Given `previousBooks` and `marketBooks` parameters, a function \f$q(\ell)\f$ that returns
         * expected quality for a given value \f$\ell\f$, and \f$\ell_{max}\f$, this numerically
         * determines the value of `l` in \f$[0, max_l]\f$ that maximizes the expected net profit
         * function:
         *
         * \f[
         *      \Pi = \beta_1 + \beta_2 (q(\ell))^D + \beta_3 firstBook + \beta_4 previousBooks +
         *      \beta_4 marketBooks - \ell
         * \f]
         *
         * using the current values of `beta`.  Note that this maximum is likely to simply be
         * \f$0\f$ or \f$\ell_{max}\f$ if \f$(q(\ell))^D\f$ is not concave.
         *
         * \param q A function, lambda or similar object that takes a double and returns the quality
         * associated with that double.  If the function returns a negative value, 0 will be
         * substituted instead.
         * \param previousBooks the previous books of the author
         * \param marketBooks the books on the market in the previous period
         * \param l_max the maximum value of `l` to consider.  This is typically the reader's
         * money (income) on hand.
         *
         * \sa eris::single_peak_search for the numerical algorithm used.
         */
        double argmaxL(
                const std::function<double(const double &)> q,
                const unsigned long &previousBooks, const unsigned long &marketBooks,
                const double &l_max
                ) const;

    private:
        unsigned int D_;

};

}}
