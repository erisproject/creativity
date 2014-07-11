#pragma once
#include <vector>
#include <unordered_map>
#include <eris/types.hpp>
#include <eris/Position.hpp>

namespace creativity {
class Book; // forward declaration
namespace state {

/** Records the various variables associated with a book.  This is basically a container class
 * with a constructor that copies the current state of a given Book. */
class BookState final {
    public:
        /** Constructs a new BookState without setting any of its values (they will be default
         * initialized).
         */
        BookState() = default;

        /// Constructs a new BookState, settings its values using the given Reader.
        BookState(const Book &r);

        /// Unique simulation ID of the reader
        eris::eris_id_t id;

        /// The simulation period this state represents.
        unsigned long t;

        /// The author of this book
        eris::eris_id_t author;

        /// Position of the reader
        eris::Position position;

        /// The quality parameter of this book, which is the mean of reader quality draws.
        double quality;

        /// True if this book was on the market this period, false if not.
        bool market;

        /// The market price of this book; will be NaN if this book was not on the market.
        double price;

        /// The revenue of this book this period.
        double revenue;

        /// The cumulative lifetime revenue of this book, up to and including the current period.
        double revenueLifetime;

        /// The number of sales of copies of this book in the current period.
        unsigned long sales;

        /// The cumulative lifetime sales of copies of this book, up to and including the current period.
        unsigned long salesLifetime;

        /** The number of copies of this book in the simulation.  This will always be at least
         * `.sales + 1` (the `+ 1` because the author has a copy of his own book), and can higher if
         * there are non-sale, i.e., piratic ways of obtaining a copy of a book.
         */
        unsigned long copies;

        /// The age of this book (in simulation periods).  0 means the book was new this period.
        unsigned long age;

        /** The number of periods this book was (or has been) on the market.  If this book is
         * currently on the market, this is simply `.age + 1`; if not on the market, this is the
         * number of periods that the book was on the market.
         */
        unsigned long lifetime;
};

}}
