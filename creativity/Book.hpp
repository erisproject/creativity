#pragma once
#include <eris/WrappedPositional.hpp>
#include <eris/Good.hpp>
#include "creativity/Reader.hpp"

namespace creativity {

// forward-declarations
class Reader;
class BookMarket;

/** Class representing a particular book.
 *
 * Lock note: an explicit lock on this class doesn't need to be held as most of it is constant.  The only
 * thing that is changeable is the number of sales, which is internally protected by sales(unsigned long).
 * Thus obtaining a lock (even a read lock) on this class is seldom necessary.
 */
class Book : public eris::Positional<eris::Good::Discrete> {
    public:
        /** Constructs a new book at the given position created by the given author.  When the book
         * is added to the simulation, a BookMarket is constructed for selling the book which
         * starts the price at `initial_price`.
         */
        Book(eris::Position p, eris::SharedMember<Reader> author, double initial_price);

        /// Returns the age of the book, in simulation periods.
        unsigned long age() const;

        /** When the Book is added to the simulation it also creates an associated BookMarket object
         * that handles sales of copies of the book.
         *
         * The creation period is stored (for use by age()), and the number of sales is initialized
         * to 0.
         */
        virtual void added() override;

        /** When the market or author is removed from the simulation, record it by clearing the
         * stored author/market fields.
         */
        void weakDepRemoved(eris::SharedMember<Member>, const eris::eris_id_t &old) override;

        /** Returns true if this Book's author is still in the simulation, false otherwise. */
        bool hasAuthor() const;

        /** Returns the author who created this Book.  If the author is no longer in the simulation,
         * this will throw an exception (via Simulation::agent()). */
        eris::SharedMember<Reader> author() const;

        /** Returns true if this book is currently on the market, false if there is no associated
         * market for this book.
         */
        bool onMarket() const;

        /** Returns the BookMarket that sells this Book.  Will throw an exception (via Simulation)
         * if no Market doesn't exist, so check onMarket() first.
         */
        eris::SharedMember<BookMarket> market() const;

        /// Returns the lifelong number of sales of this book
        unsigned long sales() const;

        /// Increases the number of sales of this book by the given amount
        void sales(unsigned long new_sales);

    private:
        unsigned long created_;
        eris::eris_id_t author_, market_;
        unsigned long copies_sold_;
        const double init_price_;
};

}

