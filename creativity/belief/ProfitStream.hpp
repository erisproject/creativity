#pragma once
#include "creativity/belief/LinearDerived.hpp"
#include "creativity/Book.hpp"
#include <eris/algorithms.hpp>
#include <unordered_set>

namespace creativity { namespace belief {

/** This class represents an author's belief about the lifetime profits of a book based on partial
 * lifetime profits.  The model is of the form:
 *
 * \f[
 *     \pi_{remaining} = \beta_0 \pi_0 + \beta_1 \pi_1 + \ldots + \beta_{K-1} \pi_{K-1} + u
 * \f]
 * where:
 * - \f$K\f$ is the number of periods of profit to consider
 * - \f$\pi_i\f$ is the profit the book earned in the period in which it had age \f$i\f$
 *
 * This model is designed to allow agents to predict the remaining lifetime profits of a book of a
 * given age \f$i \in {1, \ldots, K}\f$ by using multiple iterations of the model (with \f$K \in
 * \{1, \hdots, \overbar{K}\}\f$), so as to predict the future profitability of a book to aid in the
 * decision to keep the book on the market.
 *
 * Updating of beliefs occurs when books finish their natural life (i.e. leave the market).  At this
 * point, a book of age \f$j\f$ contributes one data point to each model with \f$K \leq j\f$: one
 * data point for each period of its life.
 *
 * \f$\beta\f$ values are not restricted.
 */
class ProfitStream : public LinearDerived<ProfitStream> {
    public:
        /** Default constructor: note that default constructed objects are not valid models.
         * \sa belief::Linear::Linear()
         */
        ProfitStream() = default;

        /** Constructs a ProfitStream with a weak prior for a model of `K` parameters.  This model
         * will have `beta` set to `0` for the first `K-1` parameters and `1` for the last
         * parameter; `s2` will be set to 1; `V` will be set to an identity matrix, and `n` will be
         * set to `1e-6`, so as to make this an extremely weak prior when used for generating a
         * posterior.
         */
        ProfitStream(unsigned int K);

        /** Constructs a ProfitStream object with the given parameter information.
         *
         * \param args a parameter pack to forward to the base Linear constructor.
         *
         * \sa Linear::Linear
         */
        template <typename ...Args>
        ProfitStream(Args &&...args) : Parent(std::forward<Args>(args)...)
        {}

        /** Given a book, this uses the profit of the first \f$K\f$ periods the book has been on the
         * market to predict the remaining cumulative lifetime profit,
         * \f$\widehat\pi_{remaining}\f$.
         *
         * \param book the book, which should have an age of at least K (otherwise the prediction
         * will be wrong as 0 will be used for revenue in \f$(age, K]\f$).
         */
        double predict(eris::SharedMember<Book> book);

        using Linear::predict;

    protected:
        /// Ensures that all beta values are non-negative
        /*virtual void verifyParameters() const override;*/

        /// Constructs a new Demand object given a Linear base object.
        virtual ProfitStream newDerived(Linear &&base) const override;

    private:
        // Initialize a ProfitStream from a Linear<K>
        explicit ProfitStream(Linear &&base) : Parent(std::move(base)) {}
};

}}
