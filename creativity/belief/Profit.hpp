#pragma once
#include "creativity/belief/Linear.hpp"
#include "creativity/Book.hpp"

namespace creativity { namespace belief {

/** This class represents an author's belief about the lifetime profitability of a work.  The model
 * is of the form:
 *
 * \f$\Pi_b = \beta_0 + \beta_1 q_b^D + \beta_2 firstBook + \beta_3 previousBooks + \beta_4 marketBooks + u\f$
 * where:
 * - \f$Pi_b\f$ is the lifetime profits of the book
 * - \f$q_b\f$ is the quality of the book.  Note that the exponentiation on this term is
 *   sign-preserving even when \f$D\f$ is even: thus for \f$q_b = -0.5, D=2\f$ the relevant quantity
 *   is \f$-0.25\f$.
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
        /** Constructs a profit model with the given parameter information.
         *
         * \param D the dimensionality of the world.
         * \param args prior arguments to forward to the base Linear constructor.
         *
         * \sa Linear::Linear
         */
        template <typename ...Args>
        Profit(unsigned int D, Args &&...args)
        : LinearBase{std::forward<Args>(args)...}, D_{D}
        {}

        /** Given a set of model parameters, this returns an expected value \f$\Pi_b\f$, the
         * lifetime profit of the book.
         *
         * \param q the quality of the book
         * \param previousBooks the number of previous books created by this book's author.  This
         * parameter also determines the `firstBook` dummy (`= 1` iff `previousBooks == 0`).
         * \param marketBooks the number of books on the market last period
         */
        double predict(double q, unsigned long previousBooks, unsigned long marketBooks) const;

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
         * \returns the double value \f$\ell\f$ that maximizes expected net profit.
         *
         * \sa eris::single_peak_search for the numerical algorithm used.
         */
        double argmaxL(
                const std::function<double(const double &)> q,
                unsigned long previousBooks, unsigned long marketBooks,
                double l_max
                ) const;

        /** Uses the current object's priors to generate a new object whose parameters are the
         * posteriors of this object after incorporating new data.
         *
         * \param y a vector of new y data
         * \param X a matrix of new X data
         */
        Profit update(const Ref<const VectorXd> &y, const Ref<const MatrixXKd> &X) const;

        /** Given a book and perceived quality, this builds an X matrix row of profit data for that
         * book.  This needs to be called after the period has advanced: typically in the
         * inter-period optimization stage.
         */
        RowVectorKd profitRow(eris::SharedMember<Book> book, double quality) const;

    private:
        // Initialize a Profit from a Linear<>
        Profit(unsigned int D, LinearBase &&base)
            : LinearBase{base}, D_{D} {}

        unsigned int D_;

};

}}
