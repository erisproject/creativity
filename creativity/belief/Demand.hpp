#pragma once
#include <eris/learning/BayesianLinearRestricted.hpp>
#include <eris/SharedMember.hpp>
#include <Eigen/Core>
#include <string>
#include <utility>

namespace creativity { class Book; }

namespace creativity { namespace belief {

/** This class represents an author's belief about the per-period demand for books.  The model is of
 * the form:
 *
 * \f$Q_b = \beta_0 + \beta_1 P_b + \beta_2 q_b + \beta_3 prevSales + \beta_4 nosales + u\f$
 * where:
 * - \f$Q_b\f$ is the quantity (i.e. copies) demanded
 * - \f$P_b\f$ is the price of a copy (which must be non-negative)
 * - \f$q_b\f$ is the quality of the work, determined when the work is created, which must be
 *   non-negative.  Readers observe either a draw based on this, or a predicted quality.
 * - \f$prevSales\f$ is the number of copies sold in previous periods
 * - \f$nosales\f$ is the number of periods the book has gone without sales
 * - \f$marketBooks\f$ is the number of books on the market in the previous period
 *
 * The following restrictions are imposed on beliefs:
 * - \f$\beta_1 \leq -0.05\f$ (demand curve is downward sloping (at least for sufficiently small p))
 * - \f$\beta_2 \geq 0\f$ (higher quality means more sales (at least for low quality))
 * - \f$\beta_4 \leq -1\f$ (each no sales period decreases quantity demanded by at least 1)
 *
 * These constraints are combined with a natural conjugate prior for the purposes of updating the
 * beliefs via Bayesian econometrics.
 *
 * Additionally, when predicting, readers know that Q can never be larger than the number of readers
 * in the world minus the copies already sold, and incorporate this limitation into predictions.
 */
class Demand : public eris::learning::BayesianLinearRestricted {
    public:
        /** Default constructor: note that unlike a default-constructed BayesianLinear model, this
         * default constructed object is a usable (but noninformative) model.
         *
         * \sa eris::learning::BayesianLinear::BayesianLinear()
         */
        Demand() : Demand(parameters()) {}

        /** Constructs a demand model with the given prior information.  This forwards arguments to
         * the BayesianLinearRestricted constructor, but also adds the model's restrictions.
         *
         * \param args prior arguments to forward to the base BayesianLinearRestricted constructor.
         *
         * \sa eris::learning::BayesianLinear::BayesianLinear
         */
        template <typename ...Args>
        explicit Demand(Args &&...args) : BayesianLinearRestricted(std::forward<Args>(args)...)
        {
            // Add restrictions:
            restrict(1) <= -0.05; // beta_price <= 0 (higher price <-> lower quantity)
            restrict(2) >= 0.0; // beta_q >= 0
            restrict(4) <= -1.0; // beta_{nosales_periods} <= 0
//            restrict(7) <= 0.0; // beta_{age} <= 0
//            restrict(9) <= 0.0; // beta_{market_size} <= 0

            // Set beta names for nicer output
            names({"const", "price", "quality", /*u8"quality²",*/ "prevSales", "noSales"});//, "I(new)", "age", "otherBooks", "marketBooks"});
        }

        /// Returns the number of parameters of this model (5)
        static unsigned int parameters() { return 5; }

        /// This model always has exactly `parameters()` parameters
        virtual unsigned int fixedModelSize() const override;

        /** Given a set of model parameters, this returns an expected value \f$Q_b\f$, the number of sales.
         *
         * \param draws the number of draws to use for prediction
         * \param P the price of a copy of the book
         * \param q the quality of the book
         * \param S prior book sales
         * \param nosales the number of time periods since the book last had a sale; if greater than
         * age, age is used (so that `sim->t() - book->lastSale()` can be passed)
         * \param age of the book
         * \param otherBooks the number of other books created by this book's author.  This parameter
         * also determines the `onlyBook` dummy (`= 1` iff `otherBooks == 0`).
         * \param lag_marketBooks the number of books on the market in the previous period
         *
         * \throws std::domain_error if `P < 0` or `q < 0`.
         */
        double predict(unsigned int draws, double P, double q, unsigned int S, unsigned int nosales, unsigned int age, unsigned int otherBooks, unsigned int lag_marketBooks);

        using BayesianLinearRestricted::predict;

        /** Given various information about a book, returns an X matrix row of data representing
         * that information.
         *
         * \param P the price of a copy of the book
         * \param q the quality of the book
         * \param S prior book sales
         * \param nosales the number of time periods since the book last had a sale; if greater than
         * age, age is used (so that `sim->t() - book->lastSale()` can be passed)
         * \param age of the book
         * \param otherBooks the number of other books created by this book's author.  This parameter
         * also determines the `onlyBook` dummy (`= 1` iff `otherBooks == 0`).
         * \param lag_marketBooks the number of books on the market in the previous period
         */
        static Eigen::RowVectorXd row(double P, double q, unsigned S, unsigned nosales, unsigned age, unsigned otherBooks, unsigned lag_marketBooks);

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
         * Note, however, that there is a limit to Q: it can never be larger than the number of
         * readers in the world who haven't read the book, in other words, the population minus the
         * copies previously sold minus the copies pirated minus 1 (an author doesn't buy his own
         * book).
         *
         * \param draws the number of draws to use for prediction
         * \param q the quality of the book
         * \param s prior book sales
         * \param z book pirated copies
         * \param n simulation readers size (including the author)
         * \param nosales the number of periods the book has gone without sales (`age` is used if
         * this is > age).
         * \param age the age of the book, in simulation periods
         * \param otherBooks the number of other books created by this book's author.  This parameter
         * also determines the `onlyBook` dummy (`= 1` iff `otherBooks == 0`).
         * \param marketBooks the number of books on the market in the last pre-prediction period
         * (including this book, if this book was on the market)
         * \param c the per-unit cost of copies.
         * \param max_price the maximum price to set (if optimum is predicted above this, reduce to
         * this)
         *
         * \returns a std::pair of double values where `.first` is the maximizing price and
         * `.second` is the quantity predicted at that price.  If optimal price is below `c` (i.e.
         * demand is 0 or negative at a price of `c`), 0 is returned for both fields.
         *
         * \throws std::domain_error if `c < 0` or `q < 0`
         */
        std::pair<double, double> argmaxP(unsigned int draws, double q, unsigned int s, unsigned int z, unsigned int n, unsigned int nosales, unsigned int age, unsigned int otherBooks, unsigned int marketBooks, double c, double max_price);

        /** Given a book and perceived quality, this builds an X matrix row of data representing
         * that book.  This method may only be called in the inter-period optimization stage, after
         * `t` has been incremented but before new books/markets have been created.  The book must
         * still be on the market: its current market price is used to build the row.
         *
         * Note that this method typically should only be used beginning in t=3: before that,
         * lag_market_books would be 0.
         */
        static Eigen::RowVectorXd bookRow(eris::SharedMember<Book> book, double quality, unsigned int lag_market_books);

        /// Returns "Demand", the name of this model.
        virtual std::string display_name() const override { return "Demand"; }
};

}}
