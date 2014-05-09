#pragma once
#include <eris/WrappedPositional.hpp>
#include <eris/Good.hpp>

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
class Book final : public eris::WrappedPositional<eris::Good::Discrete> {
    public:
        /** Constructs a new book at the given position created by the given author.  The book will
         * have the same wrapping boundaries as the author.  When the book is added to the
         * simulation, a BookMarket is constructed for selling the book which starts the price at
         * `initial_price`.  
         *
         * \param p the position of the book
         * \param author the author of the book
         * \param order the number of this book in the author's authored books list.  `order == 0`
         * indicates that this is the author's first book; `order == 12` indicates that this is the
         * author's thirteenth book.
         * \param initial_price the price of the book for the next period
         * \param quality the quality of the book
         * \param qDraw a callable object that can be called upon to produce a quality draw for a
         * given book and reader, passed as arguments.  This is often a random draw (and so
         * different readers can get different quality values).  The returned value should ideally
         * be non-negative; negative draws will be changed to 0 values.
         *
         * Note that a reference to the author is stored internally: even if the author is removed
         * from the simulation, it will still be kept from destruction by this class.
         */
        Book(eris::Position p, eris::SharedMember<Reader> author, unsigned long order, double initial_price,
                double quality, std::function<double(const Book&, const Reader&)> qDraw);

        /// Returns the age of the book, in simulation periods.
        unsigned long age() const;

        /** Returns the order of this book in the author's set of books.  A value of 0 indicates
         * that this is the author's first book, 3 indicates an author's 4th book, etc.
         */
        const unsigned long& order() const;

        /** When the Book is added to the simulation it also creates an associated BookMarket object
         * that handles sales of copies of the book.
         *
         * The creation period is stored (for use by age()), and the number of sales is initialized
         * to 0.
         */
        void added() override;

        /** When the market or author is removed from the simulation, record it by clearing the
         * stored author/market fields.
         */
        void weakDepRemoved(eris::SharedMember<Member>, const eris::eris_id_t &old) override;

        /** Returns true if this Book's author is still in the simulation, false otherwise. */
        bool livingAuthor() const;

        /** Returns the author who created this Book.  Note that the returned author might not be in
         * the simulation (check whether id() is non-zero, call hasSimulation() on the author, or
         * call this class's livingAuthor() method).
         */
        eris::SharedMember<Reader> author() const;

        /** Returns true if this book is currently on the market, false if there is no associated
         * market for this book.
         */
        bool hasMarket() const;

        /** Returns the BookMarket that sells this Book.  Will throw an exception (via Simulation)
         * if no Market doesn't exist, so check hasMarket() first.
         */
        eris::SharedMember<BookMarket> market() const;

        /// Returns the lifelong number of sales of this book
        unsigned long sales() const;

        /// Increases the number of sales of this book by the given amount
        void sales(unsigned long new_sales);

        /** Returns the actual, fixed quality value of the book.  This is meant to be used by the
         * author and the `qDraw` callable object given during construction; readers of the book
         * should instead obtain a draw using qualityDraw().
         */
        const double& quality() const;

        /** Draws a quality value for this book and the given reader.  This is intended to be called
         * once per reader *after* the reader has obtained a copy of the book; the reader should
         * then maintain the value.  (Before the reader obtains the copy, they have to use their
         * prior and/or any network information to predict a quality).
         *
         * The returned value will always be non-negative.
         */
        double qualityDraw(const Reader &reader);

    private:
        unsigned long created_, copies_sold_;
        eris::SharedMember<Reader> author_;
        const unsigned long order_;
        eris::eris_id_t market_;
        const double init_price_, quality_;
        std::function<double(const Book&, const Reader&)> quality_draw_;
};

}

