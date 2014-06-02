#pragma once
#include "creativity/belief/Linear.hpp"
#include "creativity/Book.hpp"
#include <eris/algorithms.hpp>
#include <unordered_set>

using Eigen::VectorXd;
using Eigen::MatrixXd;

namespace creativity { namespace belief {

/** This class represents an author's belief about the lifetime profits of a book based on partial
 * lifetime profits.  The model is of the form:
 *
 * \f[
 *     \pi_{remaining} = \beta_1 I_1 \pi_0 + \beta_2 I_2 (\pi_0 + \pi_1) + \beta_3 I_3 (\pi_0 +
 *     \pi_1 + \pi_2) + \ldots + \beta_{K} I_{K} (\pi_0 + \ldots + \pi_{K-1}) + u
 * \f]
 * where:
 * - \f$K\f$ is the maximum number of periods of profit to consider
 * - \f$I_i\f$ is an indicator variable that equals 1 if the book is \f$i\f$ periods old, 0 otherwise.
 * - \f$\pi_i\f$ is the profit the book earned in the period in which it had age \f$i\f$
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
        using LinearBase::Linear;

        /** Given the cumulative profit of the first \f$n\f$ periods a book is on the market, this
         * returns a predicted remaining cumulative lifetime profit, \f$\widehat\pi_{remaining}\f$.
         *
         * \param profit_curr the current cumulative profit of the book
         * \param age the age of the book
         *
         * \returns 0 if `profit_curr <= 0 || age >= K`, otherwise returns the predicted remaining
         * profit.
         */
        double predict(double profit_curr, unsigned int age) const;

        /** Uses the current object's priors to generate a new object whose parameters are the
         * posteriors of this object after incorporating new data.
         *
         * \param y a vector of new y data
         * \param X a matrix of new X data
         */
        ProfitStream update(const Ref<const VectorXd> &y, const Ref<const MatrixXKd> &X) const;

        /// Used to track which books have been added to this ProfitStream
        std::unordered_set<eris::eris_id_t> tracked;

        /** Takes a book, y/X references, and a book number, and populates the appropriate rows of y
         * and X to be passed to update().
         *
         * \param book the book to be added
         * \param y a reference to the y vector, which should be NK by 1
         * \param X a reference to the X matrix, which should be NK by K
         * \param n the index of the book; K rows beginning at nK will be populated in y and X
         */
        void populate(eris::SharedMember<Book> book, Ref<VectorXd> y, Ref<MatrixXKd> X, size_t n) const;

    protected:
        /// Ensures that all beta values are non-negative
        virtual void verifyParameters() const override;

    private:
        // Initialize a ProfitStream from a Linear<7>
        ProfitStream(LinearBase &&base) : LinearBase{base} {}
};

}}
