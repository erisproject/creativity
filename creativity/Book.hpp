#pragma once
#include <eris/WrappedPositional.hpp>
#include <eris/good/Discrete.hpp>
#include <map>
#include <functional>
#include <eris/random/truncated_normal_distribution.hpp>

namespace creativity {

// forward-declarations
class Reader;
class BookMarket;
class Creativity;

/** Class representing a particular book.
 *
 * Lock note: an explicit lock on this class doesn't need to be held as most of it is constant.  The
 * only thing that is changeable is the number of sales and revenue of those sales.
 *
 * This class automatically establishes momentary read/write locks when accessing or modifying the
 * sales and revenue information.  However, when outside locks are also in use it is recommended to
 * include the book in the encompassing lock.
 *
 * Books have three potential states: they can be on the market controlled by the author (``primary
 * market''); on the market controlled by someone else (``secondary market'') such as the
 * PublicTracker agent; or off the market.
 *
 * Without a secondary market provider such as PublicTracker, books start out in the primary market
 * and eventually move off market.  With a secondary market provider, books can start out in either
 * a primary market and then move to secondary market, or start out directly in the secondary
 * market.  See PublicTracker for details.
 */
class Book final : public eris::WrappedPositional<eris::good::Discrete> {
    public:
        /// Not default constructible
        Book() = delete;
        /** Constructs a new book at the given position created by the given author.
         *
         * \param creativity the Creativity object that owns the simulation this book belongs to
         * \param p the position of the book
         * \param author the author of the book
         * \param order the number of this book in the author's authored books list.  `order == 0`
         * indicates that this is the author's first book; `order == 12` indicates that this is the
         * author's thirteenth book.
         * \param quality the mean quality of the book; perceived quality draws are from a normal
         * with this mean and standard deviation as configured in the Creativity settings.
         *
         * Note that a reference to the author is stored internally: even if the author is removed
         * from the simulation, it will still be kept from destruction by this class.
         */
        Book(Creativity &creativity, const eris::Position &p, eris::SharedMember<Reader> author,
                unsigned int order, double quality);

        /// Returns the age of the book, in simulation periods.
        eris::eris_time_t age() const;

        /** Returns the simulation period in which the book became available (i.e. the book was
         * written at the beginning of `created()`.
         */
        const eris::eris_time_t& created() const;

        /** Returns the first simulation period in which the book was no longer on the private
         * market.  The returned value should only be used if hasPrivateMarket() returns false.
         *
         * \sa hasPrivateMarket()
         * \sa marketPeriods()
         * \sa age()
         */
        const eris::eris_time_t& leftPrivateMarket() const;

        /** Returns the simulation period in which the book was first placed on the public market or
         * 0 if the book has never been on the public market.
         */
        const eris::eris_time_t& publicMarketCreated() const;

        /** Returns the number of periods the book has been (or was) available on the private
         * market, including the current period (if the book is still on the market).  The value
         * returned by this method increases by one each period of the simulation until the book
         * exits the market, from which point on it remains fixed.
         *
         * For private-market books, this equals `age() + 1`; otherwise this equals
         * `leftPrivateMarket() - created()`.
         *
         * \sa age()
         * \sa outOfPrint()
         * \sa created()
         */
        unsigned int privateMarketPeriods() const;

        /** Returns the order of this book in the author's set of books.  A value of 0 indicates
         * that this is the author's first book, 3 indicates an author's 4th book, etc.
         */
        const unsigned int& order() const;

        /** When the Book is added to the simulation it also creates an associated BookMarket object
         * that handles sales of copies of the book.
         *
         * The creation period is stored (for use by age()), and the number of sales is initialized
         * to 0.
         */
        void added() override;

        /** Sets a new market for this book.  Throws an exception if the book already has a market,
         * or if the given market is not a market for this book.  Called automatically by
         * BookMarket; calling this externally should not be required.
         *
         * \param market the market object providing copies of this book
         */
        void setMarket(eris::SharedMember<BookMarket> market);

        /** When the market is removed from the simulation, record it by clearing the stored market.
         */
        void weakDepRemoved(eris::SharedMember<eris::Member>, eris::eris_id_t old) override;

        /** Returns true if this Book's author is still in the simulation, false otherwise. */
        bool livingAuthor() const;

        /** Returns the author who created this Book.  Note that the returned author might not be in
         * the simulation (check whether id() is non-zero, call hasSimulation() on the author, or
         * call this class's livingAuthor() method).
         */
        eris::SharedMember<Reader> author() const;

        /** Returns true if this book is currently on the market, false if there is no associated
         * market for this book.  Note that this could be the author's market, or a public market
         * (if the simulation contains a suitable public agent).
         */
        bool hasAnyMarket() const;

        /** Returns true if this book's current market exists and is a private market, that is, on
         * the market as controlled by the book's author.
         */
        bool hasPrivateMarket() const;

        /** True if the book is currently on the public market.
         */
        bool hasPublicMarket() const;

        /** Returns the BookMarket that sells this Book (which could be a PublicBookMarket if
         * hasPublicMarket() is true).  Will throw an exception (via Simulation) if no Market
         * exists, so check hasAnyMarket() first.
         */
        eris::SharedMember<BookMarket> market() const;

        /** Returns the current market price of this book.  If this book is not on the market,
         * returns a quiet NaN instead.
         */
        double price() const;

        /// Returns the lifelong number of sales (on both private and public markets) of this book
        unsigned int lifeSales() const;

        /// Returns the lifelong number of private market sales of this book
        unsigned int lifeSalesPrivate() const;

        /// Returns the lifelong number of public market sales of this book
        unsigned int lifeSalesPublic() const;

        /// Returns the number of sales of this book so far in the current period
        unsigned int currSales() const;

        /** Returns the number of sales in simulation period `t`.  Calling with `t` earlier than the
         * period this book was created is allowed: it simply returns 0.
         */
        unsigned int sales(eris::eris_time_t t) const;

        /** Returns the last simulation period `t` in which this book had at least one sale, or 0 if
         * this book has never had any sales.
         */
        eris::eris_time_t lastSale() const;

        /// Returns the lifelong number of pirated copies of this book
        unsigned int lifePirated() const;

        /// Returns the number of pirated copies of this book so far in the current period
        unsigned int currPirated() const;

        /// Returns the number of pirated copies in simulation period `t`
        unsigned int pirated(eris::eris_time_t t) const;

        /// Returns the lifelong revenue (from both private and public markets) of this book
        double lifeRevenue() const;

        /// Returns the lifelong private market revenue of this book
        double lifeRevenuePrivate() const;

        /// Returns the lifelong public revenue of this book
        double lifeRevenuePublic() const;

        /** Calculates the lifetime profit of this book from the private market, not including
         * initial creation cost, but include per-period fixed costs and unit costs.
         *
         * Note that if cost_unit or cost_market are changed, this will return what the profit would
         * have been at the current levels.
         *
         * Note that it is possible for this value to be negative (typically because the market had
         * sales insufficient to cover fixed costs).
         */
        double lifeProfitPrivate() const;

        /// Returns the revenue of this book so far in the current period
        double currRevenue() const;

        /// Returns the revenue earned by this book in simulation period `t`
        double revenue(eris::eris_time_t t) const;

        /// Returns the prize money earned by the author for this book in the current period
        double currPrize() const;

        /// Returns the prize money earned by the author for this book in simulation period `t`
        double prize(eris::eris_time_t t) const;

        /// Returns the lifetime prize money earned by the author for this book
        double lifePrize() const;

        /// Returns the lifetime number of votes for this book
        uint32_t lifeVotes() const;
        /// Returns the number of votes for this book in the current period
        uint32_t currVotes() const;
        /// Returns the number of votes for this book in the simulation period `t`
        uint32_t votes(eris::eris_time_t t) const;

        /** Returns `lifeSales() + lifePirated()` (but atomically), reflecting the number of copies
         * of the book that have been made.  This should always be one less than the total in the
         * simulation, as the author also has a copy of the book in his library which is neither a
         * sale nor a pirated copy.
         */
        unsigned int copies() const;

        /** Queries the simulation for the number of copies of this book in existence.  So long as
         * no readers are removed from the simulation, this should always be exactly one larger than
         * the value as returned by copies(), because this counts the author's original copy, while
         * copies() does not.  This method is considerably more expensive than copies() as it has to
         * examine the library of every reader in the simulation.  Call copies() instead unless you
         * have a very good reason to perform such a query.
         */
        unsigned int queryCopies() const;

        /** Increase the sales and revenue of this book for the current period.  This can safely be
         * called multiple times per period.  Both the current sales/revenue values and global
         * revenue values will be increased by the call.  If the book currently has a private
         * market, these are recorded as private sales; if the book has a public market, these are
         * public sales.  If the book has no market at all, an exception is thrown.
         *
         * \param new_sales the number of new sales to record
         * \param new_revenue the amount of new revenue to record
         * \throws std::logic_error if the book has no market
         */
        void recordSale(unsigned int new_sales, double new_revenue);

        /** Increase the number of pirated copies of this book for the current period.  This can
         * safely be called multiple times per period.  Both the current piracy count and global
         * piracy count will be increased by the call.
         */
        void recordPiracy(unsigned int new_copies);

        /** Record public money received by the author for this book.
         */
        void recordPrize(double prize);

        /** Called by a reader to cast one or more votes for this book, when under the public
         * sharing with votes policy.
         */
        void vote(uint32_t count);

        /** Returns the mean quality level determined by the author at book creation time.  This
         * value is not intended for use by readers of a book: readers rather receive a subjective
         * quality from `qualityDraw()`.  Note also that this parameter is the untruncated mean: the
         * actual quality draw will be truncated to always be non-negative, which may affect the
         * mean.
         */
        double qualityMean() const;

        /** Draws a quality value for this book.  This is intended to be called once per reader when
         * the reader obtains a copy of the book; the reader should then maintain the value.
         * (Before the reader obtains the copy, they have to use their prior to predict a quality).
         *
         * The returned value will always be non-negative.
         */
        double qualityDraw();

    private:
        Creativity &creativity_;
        eris::eris_time_t created_ = 0, left_private_market_ = 0, public_market_created_ = 0;
        uint32_t copies_private_total_ = 0, copies_pirated_total_ = 0, copies_public_total_ = 0,
                 votes_total_ = 0;
        double revenue_private_total_ = 0, revenue_public_total_ = 0, prize_total_ = 0;
        std::map<eris::eris_time_t, uint32_t> copies_sold_, copies_pirated_, votes_;
        std::map<eris::eris_time_t, double> revenue_, prize_;
        eris::SharedMember<Reader> author_;
        const unsigned int order_ = 0;
        bool market_private_ = true;
        eris::eris_id_t market_id_ = 0;
        eris::random::truncated_normal_distribution<double> quality_dist_;
};

}

