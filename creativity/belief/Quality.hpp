#pragma once
#include "creativity/belief/Linear.hpp"
#include "creativity/Book.hpp"
#include "creativity/BookMarket.hpp"
#include <eris/algorithms.hpp>
#include <eris/SharedMember.hpp>

namespace creativity {

class Book; // forward declaration

namespace belief {

/** This class represents a reader's belief about the quality of an unread book.  The model is:
 *
 * \f[
 *     q = \beta_0 + \beta_1 firstBook + \beta_2 prevBooks + \beta_3 age + \beta_4 P + \beta_5 P
 *     \times age + \beta_6 numSales + u_{i}
 * \f]
 * where:
 * - \f$firstBook\f$ is a dummy: 1 iff this was the author's first book
 * - \f$prevBooks\f$ is the number of previous books written by this author (so if this was the
 *   author's fourth book, \f$firstBook=0, prevBooks=3\f$).
 * - \f$age\f$ is the age of this book in number of periods since it was released.
 * - \f$P\f$ is the market price of the book
 * - \f$numCopies\f$ is the number of copies of the book that have been sold
 *
 * The model is updated using Bayesian econometrics as new books (and realized quality values of
 * those books) are obtained.
 */
class Quality : public Linear {
    public:
        /** Default constructor: note that default constructed objects are not valid models.
         * \sa belief::Linear::Linear()
         */
        Quality() = default;

        template <typename ...Args>
        Quality(Args &&...args) : Linear{std::forward<Args>(args)...}
        {}

        /** Given a book, this returns \f$\widehat q_b\f$, the expected quality of the book.
         *
         * \param book the book being considered.
         */
        double predict(const Book &book) const;

        /** Uses the current object's priors to generate a new object whose parameters are the
         * posteriors of this object after incorporating new data.
         *
         * \param y a vector of new y data
         * \param X a matrix of new X data
         */
        Quality update(const Eigen::Ref<const Eigen::VectorXd> &y, const Eigen::Ref<const Eigen::MatrixXd> &X) const;

        /** Given a container of books, this builds an X matrix of data representing those books.
         */
        template <class Container, typename = typename std::enable_if<std::is_same<typename Container::value_type, eris::SharedMember<Book>>::value>::type>
        Eigen::MatrixXd bookData(const Container books) {
            Eigen::MatrixXd X(books.size(), K());
            size_t i = 0;
            for (const eris::SharedMember<Book> &book : books) {
                X(i, 0) = 1;
                X(i, 1) = book->order() == 0;
                X(i, 2) = book->order();
                X(i, 3) = book->age();
                double price = book->hasMarket() ? book->market()->price() : 0.0;
                X(i, 4) = price;
                X(i, 5) = price*book->age();
                X(i, 6) = book->lifeSales();
                i++;
            }
            return X;
        }
    private:
        // Initialize a Quality from a Linear<7>
        Quality(Linear &&base) : Linear{base} {}
};

}}
