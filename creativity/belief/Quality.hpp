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
 *     q = \beta_0 + \beta_1 prevBooks + \beta_2 age + \beta_3 numCopies + u_{i}
 * \f]
 * where:
 * - \f$prevBooks\f$ is the number of previous books written by this author (so if this was the
 *   author's fourth book, \f$firstBook=0, prevBooks=3\f$).
 * - \f$age\f$ is the age of this book in number of periods since it was released.
 * - \f$numCopies\f$ is the number of copies of the book that exist in the world (sold or pirated)
 *
 * The model is updated using Bayesian econometrics as new books (and realized quality values of
 * those books) are obtained.
 */
class Quality : public Linear {
    public:
        /** Default constructor: note that unlike default constructed Linear objects, this creates a noninformative model.
         *
         * \sa belief::Linear::Linear()
         */
        Quality() : Quality(parameters()) {}

        /** Constructs a Quality object with the given parameter information.
         *
         * \param args a parameter pack to forward to the base Linear constructor.
         * 
         * \sa Linear::Linear
         */
        template <typename ...Args>
        explicit Quality(Args &&...args) : Linear(std::forward<Args>(args)...)
        {
            // Set beta names for nicer output
            names({"const", "prevBooks", "age", "numCopies"});
        }

        /// Returns the number of parameters of this model (4)
        static unsigned int parameters() { return 4; }

        /// Returns `parameters()`
        virtual unsigned int fixedModelSize() const override;

        /** Given a book, this returns \f$\widehat q_b\f$, the predicted quality of the book.
         *
         * \param book the book being considered.
         * \param draws the number of beta draws to use for prediction
         */
        double predict(const Book &book, unsigned int draws);

        using Linear::predict;

        /** Given a container of books, this builds an X matrix of data representing those books.
         */
        template <class Container, typename = typename std::enable_if<std::is_same<typename Container::value_type, eris::SharedMember<Book>>::value>::type>
        static Eigen::MatrixXd bookData(const Container books) {
            Eigen::MatrixXd X(books.size(), parameters());
            size_t i = 0;
            for (const eris::SharedMember<Book> &book : books) {
                X(i, 0) = 1;
                X(i, 1) = book->order();
                X(i, 2) = book->age();
                X(i, 3) = book->lifeSales() + book->lifePirated();
                i++;
            }
            return X;
        }

        /// Returns "Quality", the name of this model
        virtual std::string display_name() const override { return "Quality"; }

        CREATIVITY_LINEAR_DERIVED_COMMON_METHODS(Quality)
};

}}
