#pragma once
#include <eris/WrappedPositional.hpp>
#include <eris/Good.hpp>
#include <map>

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
class Book final : public eris::WrappedPositional<eris::Good::Discrete> {
    public:
        /** Constructs a new book at the given position created by the given author.
         *
         * \param creativity the Creativity object that owns the simulation this book belongs to
         * \param p the position of the book
         * \param author the author of the book
         * \param order the number of this book in the author's authored books list.  `order == 0`
         * indicates that this is the author's first book; `order == 12` indicates that this is the
         * author's thirteenth book.
         * \param quality the quality of the book
         * \param qDraw a callable object that can be called upon to produce a quality draw for a
         * given book, passed as the argument.  This is often a random draw (and so different
         * readers can get different quality values).  The returned value should ideally be
         * non-negative; negative draws will be changed to 0 values.
         *
         * Note that a reference to the author is stored internally: even if the author is removed
         * from the simulation, it will still be kept from destruction by this class.
         */
        Book(std::shared_ptr<Creativity> creativity, const eris::Position &p, eris::SharedMember<Reader> author, unsigned int order,
                double quality, std::function<double(const Book&)> qDraw);

        /// Returns the age of the book, in simulation periods.
        eris::eris_time_t age() const;

        /** Returns the simulation period in which the book became available (i.e. the book was
         * written at the beginning of `created()`.
         */
        const eris::eris_time_t& created() const;

        /** Returns the first simulation period in which the book was no longer on the market.  The
         * returned value should only be used if hasMarket() returns false.
         *
         * \sa hasMarket()
         * \sa marketPeriods()
         * \sa age()
         */
        const eris::eris_time_t& outOfPrint() const;

        /** Returns the number of periods the book has been (or was) available on the market,
         * including the current period (if the book is still on the market).  The value returned by
         * this method increases by one each period of the simulation until the book exits the
         * market, from which point on it remains fixed.
         *
         * For on-market books, this equals `age() + 1`; for off-market books, this equals
         * `outOfPrint() - created()`.
         *
         * \sa age()
         * \sa outOfPrint()
         * \sa created()
         */
        unsigned int marketPeriods() const;

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
         * \param primary true if the market is a primary (i.e. author-controlled) market, false
         * otherwise.
         */
        void setMarket(eris::SharedMember<BookMarket> market, bool primary);

        /** When the market is removed from the simulation, record it by clearing the stored market.
         */
        void weakDepRemoved(eris::SharedMember<Member>, eris::eris_id_t old) override;

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
        bool hasMarket() const;

        /** Returns true if this book is currently on the primary market, that is, on the market as
         * controlled by the book's author.  This can return false when hasMarket() returns true if
         * the book is on a secondary market, such as PublicTrackerMarket.
         */
        bool hasPrimaryMarket() const;

        /** Returns true if this book has a market, but that market is not the primary market (that
         * is, it's on a market not controlled by the author).
         */

        /** Returns the BookMarket that sells this Book.  Will throw an exception (via Simulation)
         * if no Market doesn't exist, so check hasMarket() first.
         */
        eris::SharedMember<BookMarket> market() const;

        /** Returns the current market price of this book.  If this book is not on the market,
         * returns a quiet NaN instead.
         */
        double price() const;

        /// Returns the lifelong number of sales of this book
        unsigned int lifeSales() const;

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

        /// Returns the lifelong revenue of this book
        double lifeRevenue() const;

        /// Returns the revenue of this book so far in the current period
        double currRevenue() const;

        /// Returns the revenue earned by this book in simulation period `t`
        double revenue(eris::eris_time_t t) const;

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
         * examine the library of every reader in the simulation.  Calling copies() is generally
         * preferable.
         */
        unsigned int queryCopies() const;

        /** Increase the sales and revenue of this book for the current period.  This can safely be
         * called multiple times per period.  Both the current sales/revenue values and global
         * revenue values will be increased by the call.
         */
        void recordSale(unsigned int new_sales, double new_revenue);

        /** Increase the number of pirated copies of this book for the current period.  This can
         * safely be called multiple times per period.  Both the current piracy count and global
         * piracy count will be increased by the call.
         */
        void recordPiracy(unsigned int new_copies);

        /** Returns the actual, fixed quality value of the book.  This is meant to be used by the
         * author and the `qDraw` callable object given during construction; readers of the book
         * should instead obtain a draw using qualityDraw().
         */
        const double& quality() const;

        /** Draws a quality value for this book.  This is intended to be called once per reader
         * *after* the reader has obtained a copy of the book; the reader should then maintain the
         * value.  (Before the reader obtains the copy, they have to use their prior and/or any
         * network information to predict a quality).
         *
         * The returned value will always be non-negative.
         */
        double qualityDraw();

    private:
        std::shared_ptr<Creativity> creativity_;
        eris::eris_time_t created_, out_of_print_;
        unsigned int copies_sold_total_, copies_pirated_total_;
        double revenue_total_;
        std::map<eris::eris_time_t, unsigned int> copies_sold_, copies_pirated_;
        std::map<eris::eris_time_t, double> revenue_;
        eris::SharedMember<Reader> author_;
        const unsigned int order_;
        bool market_primary_;
        eris::eris_id_t market_;
        const double quality_;
        std::function<double(const Book&)> quality_draw_;
};

}

