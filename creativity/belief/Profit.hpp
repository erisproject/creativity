#pragma once
#include <eris/learning/BayesianLinear.hpp>
#include <eris/learning/BayesianLinearRestricted.hpp>
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
 * - \f$marketBooks\f$ is the average number of books on the market in the past `creation_time`
 *   periods (so that it includes one full creation cycle).  Ideally we'd want books in the
 *   current period can't be used because it is unknown at the time the reader is making the
 *   decision of whether or not to create.
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
class Profit : public eris::learning::BayesianLinearRestricted {
    public:
        /** Construct a noninformative profit model.
         *
         * \sa eris::learning::BayesianLinear::BayesianLinear
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
         * determines the value of `l` in \f$[0, max_l]\f$ that maximizes the profit less effort.
         * Initial fixed creation cost is not included.
         *
         * \param draws the number of beta draws to use for prediction
         * \param q A function, lambda or similar object that takes a double and returns the quality
         * associated with that double.  If the function returns a negative value, 0 will be
         * substituted instead.
         * \param previousBooks the previous books of the author
         * \param marketBooks the books on the market in the previous period
         * \param l_min the minimum value of `l` to consider.  The value is typically at least the
         * effort level associated with book of `quality = MC + distancePenalty(0) +
         * sizePenalty(1)`: as long as they can observe the quality, no reader would ever buy a book
         * with lower quality.  (Note, however, that this assumes readers directly observe the
         * quality and have identical penalty functions).
         * \param l_max the maximum value of `l` to consider.  This is typically the reader's
         * money (income) on hand.
         *
         * \returns a pair of double values: the first is the double value \f$\ell\f$ that maximizes
         * expected net profit; the second is the maximized profit value *after* subtracting
         * the maximizing \f$ell\f$, but *before* subtracting fixed creation costs.
         *
         * \sa eris::single_peak_search for the numerical algorithm used.
         */
    std::pair<double, double> argmaxL(
                unsigned int draws,
                const std::function<double(const double &)> q,
                unsigned long previousBooks, unsigned long marketBooks,
                double l_min,
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
         * \param quality the mean quality of the book as determined at creation time by the book's
         * author.  Must be non-negative.
         * \param avg_market_books the average number of books on the market over the previous creation time periods
         *
         * \returns a row containing the book's information as used for this model
         */
        static Eigen::RowVectorXd profitRow(double quality, int previous_books, int avg_market_books);

        /// Returns "Profit", the name of this model
        virtual std::string display_name() const override { return "Profit"; }
};

}}
