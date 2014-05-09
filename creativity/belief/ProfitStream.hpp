#pragma once
#include "creativity/belief/Linear.hpp"
#include <eris/algorithms.hpp>

using Eigen::VectorXd;
using Eigen::MatrixXd;

namespace creativity { namespace belief {

/** This class represents an author's belief about the lifetime profits of a book based on partial
 * lifetime profits.  The model is of the form:
 *
 * \f[
 *     \pi_{remaining} = \beta_0 I_0 \pi_0 + \beta_1 I_1 (\pi_0 + \pi_1) + \beta_2 I_2 (\pi_0 +
 *     \pi_1 + \pi_2) + \ldots + \beta_{K-1} I_{t-1} (\pi_0 + \ldots + \pi_{K-1}) + u
 * \f]
 * where:
 * - \f$K\f$ is the maximum number of periods of profit to consider
 * - \f$I_i\f$ is an indicator variable that equals 1 if the book is \f$i\f$ periods old, 0 otherwise.
 * - \f$\pi_i\f$ is the profit the book earned in period \f$i\f$
 *
 * This model is designed to allow agents to predict the lifetime profits of a book of age \f$i \in {0, \ldots,
 * K-1}\f$ which has not completed its life (i.e. is still on the market); the predicted profit can
 * then be used as a data point in the lifetime profit belief.
 *
 * Updating of beliefs occurs when books finish their natural life (i.e. leave the market).  At this
 * point, each book contributes K data points: one data point for each period of its life.  For
 * example, the 4th data point for a book is:
 *
 * \f[
 *     y_4 = \pi_{remaining} = \sum_{i=4}^{K} \pi_i \\
 *     X_4 = (0, 0, 0, \sum_{i=0}^{3} \pi_i, 0, \ldots, 0)
 * \f]
 *
 * Restrictions are imposed on the model prior that all \f$\beta_i\f$ coefficients are non-negative.
 *
 * This restriction with a natural conjugate prior is used for the purposes of updating the beliefs
 * via Bayesian econometrics.
 */
class ProfitStream : public Linear<> {
    public:
        /** Constructs a demand model with the given prior information.
         *
         * \param beta_prior the prior of the mean values of beta as a column vector, which must all
         * be non-negative.  The length of this vector determines \f$K\f$, the maximum number of
         * periods to consider.
         * \param s_prior the prior of sigma (typically an estimate thereof).
         * \param V_prior the prior covariance matrix of the estimators, *without* premultiplication
         * by \f$\sigma^2\f$.  That is, for a prior from OLS, this is the matrix \f$(X^\top
         * X)^{-1}\f$, not the matrix \f$s^2(X^\top X)^{-1}\f$.  This matrix should be symmetric,
         * positive definite.
         * \param n_prior the number of data points supporting the prior (which need not be an
         * integer).
         *
         * \throws std::domain_error if any value of beta_prior is negative.
         */
        ProfitStream(
                const VectorXd &beta_prior,
                const double &s_prior,
                const MatrixXd &V_prior,
                const double &n_prior
              );

        /** Given the cumulative profit of the first \f$n\f$ periods a book is on the market, this
         * returns a predicted remaining cumulative lifetime profit, \f$\widehat\pi_{remaining}\f$.
         *
         * \param profit_curr the current cumulative profit of the book
         * \param age the age of the book
         *
         * \returns 0 if `profit_curr <= 0 || age >= K`, otherwise returns the predicted remaining
         * profit.
         */
        double predict(const double &profit_curr, const unsigned int &age) const;
};

}}
