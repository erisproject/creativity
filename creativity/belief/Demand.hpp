#pragma once
#include <creativity/belief/Linear.hpp>

namespace creativity { namespace belief {

/** This class represents an author's belief about the per-period demand for books.  The model is of
 * the form:
 *
 * \f$Q_b = \beta_1 + \beta_2 P_b^D + \beta_3 q_b^D + \beta_4 S_{b-} + \beta_5 onlyBook + \beta_6 otherBooks 
 * + \beta_7 marketBooks + u\f$
 * where:
 * - \f$Q_b\f$ is the quantity (i.e. copies) sold
 * - \f$P_b\f$ is the price of a copy (which must be non-negative)
 * - \f$q_b\f$ is the quality of the work, determined when the work is created, which must be
 *   non-negative.
 * - \f$D\f$ is the fixed dimensionality of the model (e.g. 2 for a two-dimensional world).  Both
 *   \f$q_b\f$ and \f$P_b\f$ are raised to the dimensionality because changes in either affect
 *   the radius of potential customers, with total customers being proportional to \f$r^D\f$.
 * - \f$S_{b-}\f$ is the number of copies sold in previous periods
 * - \f$onlyBook\f$ is a dummy: 1 if this is the creator's only work, 0 if the creator has other
 *   works.  (Note that this can change during a book's lifetime.)
 * - \f$otherBooks\f$ is the number of other books the author has created.  (Note that this can
 *   change during a book's lifetime.)
 * - \f$marketBooks\f$ is the number of books on the market in the previous period
 *
 * The following restrictions are imposed on beliefs:
 * - \f$\beta_2 \leq 0\f$ (higher price means fewer sales)
 * - \f$\beta_3 \geq 0\f$ (higher quality means more sales)
 * - \f$\beta_7 \leq 0\f$ (more competition means fewer sales)
 *
 * These constraints are combined with a natural conjugate prior for the purposes of updating the
 * beliefs via Bayesian econometrics.
 */
class Demand : public Linear<7> {
    public:
        /** Constructs a demand model with the given prior information.
         *
         * \param D the dimensionality of the world.
         * \param beta_prior the prior of the mean values of beta.  Must be a 7-value (column)
         * vector.  Values are in the order given in this class's description.
         * \param s_prior the prior of sigma (typically an estimate thereof).
         * \param V_prior the prior covariance matrix of the estimators, *without* premultiplication
         * by \f$\sigma^2\f$.  That is, for a prior from OLS, this is the matrix \f$(X^\top
         * X)^{-1}\f$, not the matrix \f$s^2(X^\top X)^{-1}\f$.  This matrix should be symmetric,
         * positive definite.
         * \param n_prior the number of data points supporting the prior (which need not be an
         * integer).
         */
        Demand(
                unsigned int D,
                Matrix<double, K, 1> beta_prior,
                double s_prior,
                Matrix<double, K, K> V_prior,
                double n_prior
              );

        /** Given a set of model parameters, this returns an expected value \f$Q_b\f$, the number of sales.
         *
         * \param P the price of a copy of the book
         * \param q the quality of the book
         * \param S prior book sales
         * \param otherBooks the number of other books created by this book's author.  This parameter
         * also determines the `onlyBook` dummy (`= 1` iff `otherBooks == 0`).
         * \param marketBooks the number of books on the market last period
         *
         * \throws std::domain_error if `P < 0` or `q < 0`.
         */
        double predict(const double &P, const double &q, const unsigned long &S,
                const unsigned long &otherBooks, const unsigned long &marketBooks);

        /** Given a set of model parameters (other than \f$P_b\f$) and an optional per-unit cost
         * (defaulting to 0), this returns the \f$P_b\f$ value that maximizes total profits:
         *
         * \f[
         *     P_b Q_b(P_b, \hdots) - c_b Q_b(P_b, \hdots)
         * \f]
         *
         * This optimal value is calculated numerically.
         *
         * \param q the quality of the book
         * \param S prior book sales
         * \param otherBooks the number of other books created by this book's author.  This parameter
         * also determines the `onlyBook` dummy (`= 1` iff `otherBooks == 0`).
         * \param marketBooks the number of books on the market last period
         *
         * \sa eris::single_peak_search for the numerical algorithm used.
         */
        double argmax(const double &q, const unsigned long &S, const unsigned long &otherBooks, const unsigned long &marketBooks,
                const double &c = 0.0);
        // FIXME: implement this

    private:
        unsigned int D_;

};

}}
