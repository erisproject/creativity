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
        using LinearBase::Linear;

        /** Given a book, this returns \f$\widehat q_b\f$, the expected quality of the book.
         *
         * \param book the book being considered.
         */
        double predict(const Book &book) const;

        /** Uses the current object's priors to generate a new object whose parameters are the
         * posteriors of this object after incorporating new data.
         *
         * \param y a vector of new y data
         * \param X a matrix of new X data
         */
        Quality update(const Ref<const VectorXd> &y, const Ref<const MatrixXKd> &X) const;

    private:
        // Initialize a Quality from a Linear<7>
        Quality(LinearBase &&base) : LinearBase{base} {}
};

}}
