#pragma once
#include "creativity/belief/LinearRestricted.hpp"
#include "creativity/belief/LinearDerived.hpp"
#include <eris/algorithms.hpp>
#include <eris/SharedMember.hpp>

#include <eris/debug.hpp>

namespace creativity {
class Book;
namespace belief {

/** This class represents an author's belief about the per-period demand for books.  The model is of
 * the form:
 *
 * \f$Q_b = \beta_0 + \beta_1 P_b + \beta_2 q_b + \beta_3 q_b^2 + \beta_4 S_{b-} +
 * \beta_5 age + \beta_6 onlyBook + \beta_7 otherBooks + \beta_8 marketBooks + u\f$
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
 * - \f$\beta_2 \geq 0\f$ (higher quality means more sales (at least for low quality))
 * - \f$\beta_3 \leq 0\f$ (quality increase effect is concave)
 * - \f$\beta_8 \leq 0\f$ (more competition means lower (individual) demand)
 *
 * These constraints are combined with a natural conjugate prior for the purposes of updating the
 * beliefs via Bayesian econometrics.
 */
class Demand : public LinearDerived<Demand, LinearRestricted> {
    public:
        /** Default constructor: note that default constructed objects are not valid models.
         * \sa belief::Linear::Linear()
         */
        Demand() = default;

        /** Constructs a noninformative demand model.
         *
         * \param D the dimensionality of the world.
         *
         * \sa Linear::Linear
         */
        explicit Demand(unsigned int D) : Demand(D, parameters()) {}

        /** Constructs a demand model with the given prior information.
         *
         * \param D the dimensionality of the world.
         * \param args prior arguments to forward to the base Linear constructor.
         *
         * \sa Linear::Linear
         */
        template <typename ...Args>
        explicit Demand(unsigned int D, Args &&...args)
        : Parent(std::forward<Args>(args)...), D_{D}
        {
            // Add restrictions:
            restrict(1) <= 0.0; // beta_price <= 0 (higher price <-> lower quantity)
            restrict(2) >= 0.0; // beta_q >= 0
            restrict(3) <= 0.0; // beta_{q^2} <= 0
            restrict(8) <= 0.0; // more competition <-> lower demand
        }

        /// Returns the number of parameters of this model
        static unsigned int parameters() { return 9; } // NB: when changing this, change following doc, too!

        /// This model always has exactly 9 parameters
        virtual unsigned int fixedModelSize() const override;

        /** Given a set of model parameters, this returns an expected value \f$Q_b\f$, the number of sales.
         *
         * \param P the price of a copy of the book
         * \param q the quality of the book
         * \param S prior book sales
         * \param age of the book
         * \param otherBooks the number of other books created by this book's author.  This parameter
         * also determines the `onlyBook` dummy (`= 1` iff `otherBooks == 0`).
         * \param marketBooks the number of books on the market last period
         *
         * \throws std::domain_error if `P < 0` or `q < 0`.
         */
        double predict(double P, double q, unsigned long S, unsigned long age, unsigned long otherBooks, unsigned long marketBooks);

        using LinearRestricted::predict;

        /** Given a set of model parameters (other than \f$P_b\f$) and a per-unit cost this returns
         * the \f$P_b\f$ value that maximizes expected total profits:
         *
         * \f[
         *     P_b Q_b(P_b, \ldots) - c_b Q_b(P_b, \ldots)
         * \f]
         * which, due to the linearity of P_b in Q_b, has a maximum at:
         * \f[
         *     P_b^* = \frac{X\beta_{-1}}{-2\beta_1} + \frac{c_b}{2}
         * \f]
         * where \f$X\f$ is all of the model terms that do not depend on \f$P_b\f$ and
         * \f$\beta_{-1}\f$ are the associated coefficients.
         *
         * \param q the quality of the book
         * \param S prior book sales
         * \param age of the book
         * \param otherBooks the number of other books created by this book's author.  This parameter
         * also determines the `onlyBook` dummy (`= 1` iff `otherBooks == 0`).
         * \param marketBooks the number of books on the market in the last pre-prediction period
         * (including this book, if this book was on the market)
         * \param c the per-unit cost of copies.
         *
         * \returns a std::pair of double values where `.first` is the maximizing price and
         * `.second` is the maximum at that price.  If optimal profits are negative (i.e. there is
         * no demand even at minimum price `c`), 0 is returned for both fields.
         *
         * \throws std::domain_error if `c < 0` or `q < 0`
         */
        std::pair<double, double> argmaxP(double q, unsigned long S, unsigned long age, unsigned long otherBooks, unsigned long marketBooks, double c);

        /** Given a book and perceived quality, this builds an X matrix row of data representing
         * that book.  This method may only be called in the inter-period optimization stage, after
         * `t` has been incremented but before new books/markets have been created.  The book must
         * still be on the market: its current market price is used to build the row.
         */
        Eigen::RowVectorXd bookRow(eris::SharedMember<Book> book, double quality) const;

        /// Constructs a new Demand object given a Linear base object.
        virtual Demand newDerived(Linear &&base) const override;

    private:
        unsigned int D_;

        // Initialize a Demand from a Linear object of the correct size
        Demand(unsigned int D, Linear &&base) : Parent(std::move(base)), D_{D} {}
};

}}
