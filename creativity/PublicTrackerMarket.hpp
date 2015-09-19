#pragma once
#include <eris/Market.hpp>
#include <eris/Optimize.hpp>
#include "creativity/BookMarket.hpp"

namespace creativity {

class Book;
class Creativity;

/** This class extends BookMarket by always pricing the book at marginal cost (and thus never has
 * any profits).  It is created by PublicTracker for books that are taken off (or are never placed
 * on) the market by authors.
 */
class PublicTrackerMarket : public BookMarket, public virtual eris::interopt::Advance {
    public:
        /** Constructs a PublicTrackerMarket that sells copies of the given Book at marginal cost of
         * the lower of `unit_cost` or `piracy_cost`.
         */
        PublicTrackerMarket(std::shared_ptr<Creativity> creativity, eris::SharedMember<Book> b);
        /** Sets the initial price and registers this as the book's market when added to the
         * simulation.
         */
        virtual void added() override;
        /** Updates the market price of the book to the minimum of unit_cost or piracy_cost (in case
         * those changed for the upcoming period).
         */
        virtual void interAdvance() override;
        /** Overrides BookMarket finish to just discard proceeds instead of transferring them to the
         * author: since price is at marginal cost, net proceeds will come out to 0 anyway, so there
         * is never any profit to redistribute.  (Technically, the proceeds should be transferred to
         * the PublicTracker, which should also subtract costs, but we just skip that step by just
         * clearing any proceeds).
         */
        virtual void intraFinish() override;

        /** Returns true: this is a public market. */
        virtual bool isPublic() const override { return true; }

    protected:
        /// Resets price to the smaller of cost_unit and cost_piracy.
        void updatePrice();
};

}
