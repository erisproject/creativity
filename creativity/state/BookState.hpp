#pragma once
#include <eris/Position.hpp>
#include <eris/SharedMember.hpp>
#include <cmath>

namespace creativity {
class Book; // forward declaration
namespace state {

/** Records the various variables associated with a book.  This is basically a container class
 * with a constructor that copies the current state of a given Book. */
class BookState final {
    public:
        /// Constructs a new BookState, setting its values using the given Book member.
        BookState(const eris::SharedMember<Book> &b);

        /// Constructs a new BookState, setting its value using the given Book reference.
        BookState(const Book &b);

        /** Constructs a new blank BookState for a book with a position of the given number of
         * dimensions.  All values will be default initialized.  (The number of dimensions is needed
         * for Position initialization).
         */
        explicit BookState(const unsigned int dimensions);

        /// Unique simulation ID of the book
        eris::eris_id_t id;

        /// The author of this book
        eris::eris_id_t author;

        /// Position of the book
        eris::Position position;

        /// The quality parameter of this book, which is the mean of reader quality draws.
        double quality;

        /** Returns true if this book was on either the private or public market this period, false
         * if not.  The is determined directly from the price: a NaN price indicates a non-market
         * status.
         */
        bool market_any() const { return not std::isnan(price); }

        /** True if this book is on the private market, false if on the public market or not on the
         * market at all.
         */
        bool market_private;

        /// Returns true if this book is on the public market this period.
        bool market_public() const { return market_any() and not market_private; }

        /// The market price of this book; will be NaN if this book was not on the market.
        double price;

        /// The private market revenue of this book this period.
        double revenue;

        /// The cumulative lifetime private market revenue of this book, up to and including the current period.
        double revenue_lifetime;

        /// The public/prize money received by the book's author this period for this book.
        double prize;

        /// The lifetime public/prize money received by the book's author for this book.
        double prize_lifetime;

        /** The number of sales of copies of this book in the current period.  The sales are private
         * sales if `.market_private` is true and public sales if `.market_public()` is true.
         */
        uint32_t sales;

        /// The cumulative lifetime private market sales of this book, up to and including the current period.
        uint32_t sales_lifetime_private;

        /// The cumulative lifetime public market sales of this book, up to and including the current period.
        uint32_t sales_lifetime_public;

        /** Returns the cumulative lifetime sales of this book, up to and including the current
         * period, of both private and public markets.
         */
        uint32_t sales_lifetime() const { return sales_lifetime_private + sales_lifetime_public; }

        /// The number of pirated copies of this book in the current period.
        uint32_t pirated;

        /// The cumulative lifetime pirated copies of the book, up to and including the current period.
        uint32_t pirated_lifetime;

        /** Returns the number of new copies created in the current period.  This is simply
         * `sales + pirated`.
         */
        uint32_t copies() const { return pirated + sales; }

        /** Returns the cumulative lifetime copies of the book, up to and including the current
         * period.  This is simply `sales_lifetime_private + sales_lifetime_public + pirated_lifetime`.
         */
        uint32_t copies_lifetime() const { return pirated_lifetime + sales_lifetime_private + sales_lifetime_public; }

        /// The period in which this book ws created.  Directly related to age.
        eris::eris_time_t created;

        /** The number of periods this book was (or has been) on the private market.  If this book
         * is currently on the private market, this is simply the age plus 1; if not on the market,
         * this is the number of periods that the book was on the market.
         */
        uint32_t lifetime_private;

        /// The number of votes the book received this period.
        uint32_t votes;

        /// The lifetime number of votes the book has received.
        uint32_t votes_lifetime;
};

}}
