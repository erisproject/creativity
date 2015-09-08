#pragma once
#include <eris/belief/BayesianLinear.hpp>
#include <eris/belief/BayesianLinearRestricted.hpp>
#include <Eigen/Core>
#include <algorithm>
#include <functional>
#include <string>

namespace creativity { namespace belief {

/** This class represents an author's belief about the lifetime profitability of a work.  The model
 * is:
 *
 * \f$\pi_b = \beta_0 + \beta_1 q_b + \beta_2 q_b^2 + \beta_3 marketBooks + u\f$
 *
 * Specifically, the fields above are:
 * - \f$\pi_b\f$ is the lifetime profits of the book *without* considering the initial creation effort
 * - \f$q_b\f$ is the (non-negative) quality of the book.
 * - \f$marketBooks\f$ is the number of books on the market in the previous period (books in the
 *   current period can't be used because it is unknown at the time the reader is making the
 *   decision of whether or not to create).
 *
 * Although a \f$\beta_3 \leq 0\f$ restriction would make some intuitive sense (more competition
 * means lower profit), it isn't imposed because doing so induces volatile creation cycles: readers decide
 * to create in low-book periods, which results in an immediate surge in on-market books, which
 * induces readers to not create in the next period, and so on.  Not imposing it still often has it
 * be significant for many readers, but seems to avoid the strong behavioural cycles.
 *
 * The model uses a natural conjugate prior for the purposes of updating the beliefs via Bayesian
 * econometrics.
 */
class Profit : public eris::belief::BayesianLinearRestricted {
    public:
        /** Construct a noninformative profit model.
         *
         * \sa eris::belief::BayesianLinear::BayesianLinear
         */
        Profit() : Profit(parameters()) {}

        /** Constructs a profit model with the given parameter information.
         *
         * \param args prior arguments to forward to the base BayesianLinear constructor.
         *
         * \sa BayesianLinear::BayesianLinear
         */
        template <typename ...Args>
        explicit Profit(Args &&...args) : BayesianLinearRestricted(std::forward<Args>(args)...)
        {
            // Add restrictions:
            //restrict(1) >= 0; // beta_q >= 0 (higher quality <-> higher profits, at least for low quality)
            //restrict(2) <= 0; // beta_{q^2} <= 0 (quality effect is concave)
            //restrict(3) <= 0; // beta_{marketbooks} <= 0 (more competition <-> lower profit)

            // Set beta names for nicer output
            names({"const", "quality", u8"quality²", "marketBooks"});
        }

        /// Returns the number of parameters of this model (4)
        static unsigned int parameters() { return 4; }

        /// Returns `parameters()`
        virtual unsigned int fixedModelSize() const override;

        /** Given a set of model parameters, this returns an expected value \f$\Pi_b\f$, the
         * lifetime profit of the book.
         *
         * \param draws the number of draws to use for prediction
         * \param q the quality of the book
         * \param previousBooks the number of previous books created by this book's author.  This
         * parameter also determines the `firstBook` dummy (`= 1` iff `previousBooks == 0`).
         * \param marketBooks the number of books on the market last period
         */
        double predict(unsigned int draws, double q, unsigned long previousBooks, unsigned long marketBooks);

        using BayesianLinearRestricted::predict;

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
         * \param draws the number of beta draws to use for prediction
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
                unsigned int draws,
                const std::function<double(const double &)> q,
                unsigned long previousBooks, unsigned long marketBooks,
                double l_max
                );

        /** Builds an X matrix row of profit data for a book.  This needs to be called after the
         * period has advanced: typically in the inter-period optimization stage.
         *
         * This method should only be called in periods \f$t=3\f$ and later (because otherwise the
         * lagged value is not useful: the first books are created in \f$t=1\f$, and so the lagged
         * value will be 0 at the beginning of \f$t=2\f$, thus need to wait until the beginning of
         * \f$t=3\f$).
         *
         * \param previous_books the number of previous books written by the same author
         * \param quality the quality of the book as perceived by the reader this belief is for
         * \param lag_market_books the number of books that was in the market in the period just
         * before the just-ended period.
         *
         * \returns a row containing the book's information as used for this model
         */
        static Eigen::RowVectorXd profitRow(double quality, int previous_books, int lag_market_books);

        /// Returns "Profit", the name of this model
        virtual std::string display_name() const override { return "Profit"; }

        ERIS_BELIEF_BAYESIANLINEAR_DERIVED_COMMON_METHODS(Profit)
};

}}
