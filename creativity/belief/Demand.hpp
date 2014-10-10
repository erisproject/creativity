#pragma once
#include "creativity/belief/LinearRestricted.hpp"
#include <eris/algorithms.hpp>
#include <eris/SharedMember.hpp>

#include <eris/debug.hpp>

namespace creativity {
class Book;
namespace belief {

/** This class represents an author's belief about the per-period demand for books.  The model is of
 * the form:
 *
 * \f$Q_b = \beta_0 + \beta_1 P_b + \beta_2 P_b^2 + \beta_3 q_b + \beta_4 q_b^2 + \beta_5 S_{b-} +
 * \beta_6 age + \beta_7 onlyBook + \beta_8 otherBooks + \beta_9 marketBooks + u\f$
 * where:
 * - \f$Q_b\f$ is the quantity (i.e. copies) sold
 * - \f$P_b\f$ is the price of a copy (which must be non-negative)
 * - \f$q_b\f$ is the quality of the work, determined when the work is created, which must be
 *   non-negative.
 * - \f$S_{b-}\f$ is the number of copies sold in previous periods
 * - \f$age\f$ is the age of the book, in simulation periods, starting from 0.
 * - \f$onlyBook\f$ is a dummy: 1 if this is the creator's only work, 0 if the creator has other
 *   works.  (Note that this can change during a book's lifetime.)
 * - \f$otherBooks\f$ is the number of other books the author has created.  (Note that this can
 *   change during a book's lifetime.)
 * - \f$marketBooks\f$ is the number of books on the market in the previous period
 *
 * The following restrictions are imposed on beliefs:
 * - \f$\beta_1 \leq 0\f$ (demand curve is downward sloping (at least for sufficiently small p))
 * - \f$\beta_3 \geq 0\f$ (higher quality means more sales (at least for low quality))
 * - \f$\beta_4 \leq 0\f$ (quality increase effect is concave)
 * - \f$\beta_9 \leq 0\f$ (more competition means lower (individual) demand)
 *
 * These constraints are combined with a natural conjugate prior for the purposes of updating the
 * beliefs via Bayesian econometrics.
 */
class Demand : public LinearRestricted {
    public:
        /** Default constructor: note that default constructed objects are not valid models.
         * \sa belief::Linear::Linear()
         */
        Demand() = default;

        /** Constructs a demand model with the given prior information.
         *
         * \param D the dimensionality of the world.
         * \param args prior arguments to forward to the base Linear constructor.
         *
         * \sa Linear::Linear
         */
        template <typename ...Args>
        Demand(unsigned int D, Args &&...args)
        : LinearRestricted(std::forward<Args>(args)...), D_{D}
        {
            // Add restrictions:
            upperBounds()[1] = 0.0; // beta_price <= 0 (higher price <-> lower quantity)
            lowerBounds()[3] = 0.0; // beta_q >= 0
            upperBounds()[4] = 0.0; // beta_{q^2} <= 0
            upperBounds()[9] = 0.0; // more competition <-> lower demand
        }

        /// Returns the number of parameters of this model
        static unsigned int parameters() { return 10; }

        /// This model always has exactly 8 parameters
        virtual unsigned int fixedModelSize() const override;

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
        double predict(double P, double q, unsigned long S, unsigned long otherBooks, unsigned long marketBooks);

        using LinearRestricted::predict;

        /** Given a set of model parameters (other than \f$P_b\f$) and an optional per-unit cost
         * (defaulting to 0), this returns the \f$P_b\f$ value that maximizes total profits:
         *
         * \f[
         *     P_b Q_b(P_b, \ldots) - c_b Q_b(P_b, \ldots)
         * \f]
         *
         * If \f$c_b = 0\f$ this value is calculated analytically as:
         * \f[
         *     P_b = \sqrt[D]{\frac{X\gamma}{-\beta_1(1+D)}}
         * \f]
         * where \f$X\gamma\f$ is the set of values and parameters in the model other than \f$P\f$.
         * When \f$X\gamma <= 0\f$, \f$P_b = 0\f$ is returned instead.  Note that \f$\beta_1\f$ is
         * negative by assumption (and prior) so that the equation above has the same sign as
         * \f$X\gamma\f$.
         *
         * When \f$c_b > 0\f$, the profit equation above equals:
         * \f[
         *     (P - c)X\gamma + \beta_1(P-c)P^D
         * \f]
         * where, as above, \f$X\gamma\f$ is the value of all non-price independent variables and
         * coefficients.  When \f$X\gamma \leq 0\f$, expected quantity is 0 and so the minimum
         * price, \f$P = c\f$, is returned.  Otherwise the optimum is an equation of the form
         * \f$a + bP^{D-1} + cP^D = 0\f$, which has no simple analytical solution for an arbitrary
         * \f$D\f$ dimensionality value.
         *
         * Instead the simple eris::single_peak_search numerical algorithm is used.  Note that in
         * the profit equation above, the first term (\f$(P-c)X\gamma\f$) is positive and the second
         * term (\f$\beta_1(P-c)P^D\f$) is negative since \f$\beta_1 < 0\f$ (by assumption and
         * prior).  Thus the profit equation is clearly single-peaked in \f$P\f$, and so
         * eris::single_peak_search will work without issue.
         *
         * \param q the quality of the book
         * \param S prior book sales
         * \param otherBooks the number of other books created by this book's author.  This parameter
         * also determines the `onlyBook` dummy (`= 1` iff `otherBooks == 0`).
         * \param marketBooks the number of books on the market in the last pre-prediction period
         * (including this book, if this book was on the market)
         * \param c the per-unit cost of copies.  Optional: defaults to 0.
         *
         * \returns a std::pair of double values where `.first` is the maximizing price and
         * `.second` is the maximum at that price.  The returned price will always be greater than
         * or equal to `c` (which defaults to 0).
         *
         * \throws std::domain_error if `c < 0` or `q < 0`
         *
         * \sa argmaxP_MAX: the maximum P that will be considered when using the numerical algorithm
         * (i.e. when `c > 0`).
         * \sa eris::single_peak_search for the numerical algorithm used.
         */
        std::pair<double, double> argmaxP(double q, unsigned long S, unsigned long otherBooks, unsigned long marketBooks, double c = 0.0);

        /** The maximum P that will be considered for argmaxP() when called with `c > 0`.  Only
         * needs to be adjusted if the optimal value of P could potential exceed the default of 10,000.
         */
        double argmaxP_MAX = 10000;

        /** Given a book and perceived quality, this builds an X matrix row of data representing
         * that book.  This needs to be called after the period has advanced: typically in the
         * inter-period optimization stage.
         */
        Eigen::RowVectorXd bookRow(eris::SharedMember<Book> book, double quality) const;

        /** Uses the current object's priors to generate a new object whose parameters are the
         * posteriors of this object after incorporating new data.
         *
         * \param y a vector of new y data
         * \param X a matrix of new X data
         */
        Demand update(const Eigen::Ref<const Eigen::VectorXd> &y, const Eigen::Ref<const Eigen::MatrixXd> &X) const;

    private:
        unsigned int D_;

        // Initialize a Demand from a Linear object of the correct size
        Demand(unsigned int D, Linear &&base) : LinearRestricted(std::move(base)), D_{D} {}
};

}}
