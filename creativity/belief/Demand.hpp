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
 *
 * \f$\beta_8 \leq 0\f$ seems an intuitive restriction, but isn't imposed because we actually use
 * lagged market size, not current market size, because lagged market size is all the reader has
 * when doing demand prediction.  Positive values aren't entirely possible there, in particular if
 * the market exhibits cyclical behaviour.
 *
 * These constraints are combined with a natural conjugate prior for the purposes of updating the
 * beliefs via Bayesian econometrics.
 */
class Demand : public LinearRestricted {
    public:
        /** Default constructor: note that unlike a default-constructed Linear model, this default
         * constructed object is a usable (but noninformative) model.
         *
         * \sa belief::Linear::Linear()
         */
        Demand() : Demand(parameters()) {}

        /** Constructs a demand model with the given prior information.  This forwards arguments to
         * the LinearRestricted constructor, but also adds the model's restrictions.
         *
         * \param args prior arguments to forward to the base LinearRestricted constructor.
         *
         * \sa Linear::Linear
         */
        template <typename ...Args>
        explicit Demand(Args &&...args) : LinearRestricted(std::forward<Args>(args)...)
        {
            // Add restrictions:
            restrict(1) <= 0.0; // beta_price <= 0 (higher price <-> lower quantity)
            restrict(2) >= 0.0; // beta_q >= 0
            restrict(3) <= 0.0; // beta_{q^2} <= 0
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
         * \param lag_marketBooks the number of books on the market in the previous period
         *
         * \throws std::domain_error if `P < 0` or `q < 0`.
         */
        double predict(double P, double q, unsigned int S, unsigned int age, unsigned int otherBooks, unsigned int lag_marketBooks);

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
        std::pair<double, double> argmaxP(double q, unsigned int S, unsigned int age, unsigned int otherBooks, unsigned int marketBooks, double c);

        /** Given a book and perceived quality, this builds an X matrix row of data representing
         * that book.  This method may only be called in the inter-period optimization stage, after
         * `t` has been incremented but before new books/markets have been created.  The book must
         * still be on the market: its current market price is used to build the row.
         *
         * Note that this method typically should only be used beginning in t=3: before that,
         * lag_market_books would be 0.
         */
        Eigen::RowVectorXd bookRow(eris::SharedMember<Book> book, double quality, unsigned int lag_market_books) const;

        /// Returns "Demand", the name of this model.
        virtual std::string display_name() const override { return "Demand"; }

        CREATIVITY_LINEAR_DERIVED_COMMON_METHODS(Demand)
};

}}
