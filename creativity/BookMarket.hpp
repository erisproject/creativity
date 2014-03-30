#pragma once
#include <eris/Market.hpp>
#include "creativity/Book.hpp"
#include "creativity/common.hpp"

namespace creativity {

/** A BookMarket is very simple: it has an exogenously determined price and can supply infinite
 * copies of a Book at that price.
 */
class BookMarket : public eris::Market, public virtual eris::interopt::Advance {
    public:
        /** Constructs a BookMarket that sells copies of the given Book with an initial price of
         * MONEY times p.  Books are created as needed at zero cost.
         */
        BookMarket(const Book &b, const double &price);

        /** Returns price info.  Since price is constant, and there is no quantity limits,
         * this price into is simple: it's always feasible, and total is just quantity times price,
         * truncated to an integer (since books are discrete).
         */
        price_info price(double q) const override;

        /** Returns quantity info for a given payment.  Since price is constant, this is simple:
         * .quantity is the floor of p divided by current price, .constrained is always false,
         * .spent equals p*.quantity, and .unspent equals whatever amount is leftover which can't
         * buy another whole book.
         */
        quantity_info quantity(double p) const override;

        /** Sets the price of copies of the book on the market for future sales.  The value is a
         * multiple of creativity::MONEY.
         *
         * This may not be called during an intra-optimization phase.
         */
        void setPrice(double p);

        /** Returns the price of a single copy of the book.  This will be the same value as any of
         * `m.price(1).total`, `m.price(1).marginalFirst`, `m.price(1).marginal`.
         *
         * This is guaranteed not to change during an intra-period optimization phase (and so a read
         * lock on the market to obtain market price is not required).
         */
        const double& price();

        /** Reserves q units, paying at most p_max for them.  Note that if q is not an integer, this
         * will actually purchase floor(q) units.
         */
        virtual Reservation reserve(
                eris::SharedMember<eris::AssetAgent> agent,
                double q,
                double p_max = std::numeric_limits<double>::infinity()) override;

        /** Returns the Book sold by this market. */
        eris::SharedMember<Book> book();

        /** Transfers all sales generated in the previous period to the author.
         */
        void interAdvance() override;

    protected:
        /** We override reservation completion to transfer the payment to the author and put the
         * book into the b_ bundle.
         *
         * \todo this is simple: eventually it would probably be better to have a firm (perhaps
         * created when the Book is created?) that accepts payments, mitigating this override.
         */
        virtual void buy_(Reservation_ &res) override;

    private:
        eris::eris_id_t book_;
        double price_;
        // The sales of the book (these get transferred in the interApply phase)
        eris::Bundle proceeds_;
};

}
