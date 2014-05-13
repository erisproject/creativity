#pragma once
#include "creativity/belief/Linear.hpp"
#include <eris/algorithms.hpp>

namespace creativity {

class Book; // forward declaration

namespace belief {

/** This class represents a reader's belief about the quality of an unread book.  The model is:
 *
 * \f[
 *     q = \beta_0 + \beta_1 firstBook + \beta_2 prevBooks + \beta_3 age + \beta_4 P + \beta_5 P
 *     \times age + \beta_6 numSales + u_{i}
 * \f]
 * where:
 * - \f$firstBook\f$ is a dummy: 1 iff this was the author's first book
 * - \f$prevBooks\f$ is the number of previous books written by this author (so if this was the
 *   author's fourth book, \f$firstBook=0, prevBooks=3\f$).
 * - \f$age\f$ is the age of this book in number of periods since it was released.
 * - \f$P\f$ is the market price of the book
 * - \f$numCopies\f$ is the number of copies of the book that have been sold
 *
 * The model is updated using Bayesian econometrics as new books (and realized quality values of
 * those books) are obtained.
 */
class Quality : public Linear<7> {
    public:
        /** Constructs a quality model with the given prior information.
         *
         * \param beta_prior the prior of the mean values of beta.  Must be a 7-value (column)
         * vector.  Values are in the order given in this class's description.
         * \param s2_prior the prior of \f$s^2\f$ (typically an estimate thereof).
         * \param V_prior the prior covariance matrix of the estimators, *without* premultiplication
         * by \f$\sigma^2\f$.  That is, for a prior from OLS, this is the matrix \f$(X^\top
         * X)^{-1}\f$, not the matrix \f$s^2(X^\top X)^{-1}\f$.  This matrix should be symmetric,
         * positive definite.
         * \param n_prior the number of data points supporting the prior (which need not be an
         * integer).
         */
        Quality(
                const VectorKd &beta_prior,
                const double &s2_prior,
                const MatrixKd &V_prior,
                const double &n_prior
              );

        /** Given a book, this returns \f$\widehat q_b\f$, the expected quality of the book.
         *
         * \param book the book being considered.
         */
        double predict(const Book &book) const;

};

}}
