#pragma once
#include <eris/types.hpp>

namespace creativity {

/** Container class storing the data associated with a reader's copy of a book. */
class BookCopy final {
    public:
        /// The status a book copy can have
        enum class Status {
            wrote, ///< This is the author's own copy of his book
            purchased, ///< The reader bought this book on the market
            pirated ///< The reader pirated this book
        };

        /// Creates a BookCopy
        BookCopy(double quality, Status status, eris::eris_time_t acquired)
            : quality(quality), status(status), acquired(acquired) {}

        const double quality; ///< The reader's perceived quality of the book
        const Status status; ///< The status of the book (wrote, purchased, or pirated)
        const eris::eris_time_t acquired; ///< The simulation period when this copy was acquired

        /// Alias for `status == BookCopy::Status::wrote`
        bool wrote() const { return status == Status::wrote; }

        /// Alias for `status == BookCopy::Status::pirated`
        bool pirated() const { return status == Status::pirated; }

        /// Alias for `status == BookCopy::Status::wrote`
        bool purchased() const { return status == Status::purchased; }

};

}
