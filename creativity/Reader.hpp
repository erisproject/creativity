#pragma once
#include <eris/WrappedPositional.hpp>
#include <eris/agent/AssetAgent.hpp>
#include <eris/Market.hpp>
#include <string>
#include <vector>
#include <unordered_set>

namespace creativity {

class Book; // forward declaration
class BookMarket;

/** A Reader is an agent with a position whose utility is determined by books and an outside option.
 * In particular, his utility is quasilinear of the form:
 *
 * \f$u(m, x_0, x_1, \hdots, x_{n-1}) = m + \sum_{i=1}^{n} f(b_i) - p(n)\f$
 *
 * where:
 *     - \f$r\f$ is the Reader's position
 *     - \f$b_1, \hdots, b_n\f$ are the positions of the \f$n\f$ books the agent buys this period
 *     - \f$f(b)\f$ is a function that maps a book into a utility value using a decreasing function
 *     of the distance between the reader's position and the book's position.  It returns 0 for any
 *     books that the reader already owns.
 *     - \f$p(n)\f$ is a penalty function that increases in the number of books read
 *     - \f$m\f$ is the quantity of non-book spending the reader engages in which delivers a
 *     constant marginal utility of 1.
 *
 * The reader optimizes by looking at all Books currently available for sale, assessing their
 * quality, considering their price, then deciding which books to buy.  Books are only purchased
 * once, and at most one copy is purchased by any reader.
 */
class Reader : public eris::WrappedPositional<eris::agent::AssetAgent>,
    public virtual eris::interopt::Apply,
    public virtual eris::interopt::Advance,
    public virtual eris::intraopt::OptApplyReset,
    public virtual eris::intraopt::Initialize
{
    public:
        /// Constructor takes the reader position and the (wrapping) positional boundaries
        Reader(const eris::Position &pos, const eris::Position &b1, const eris::Position &b2);

        /** Takes a money value and a container of SharedMember<Book> objects and returns the
         * utility for that set of books.  This is calculated as the money plus the sum of utilities
         * from each book (from uBook()) minus the penalty incurred by buying the given number of
         * books.
         *
         * \param the money leftover after buying the provided set of books
         * \param books an iterable container of SharedMember<Book> objects.
         */
        template <typename Container>
        typename std::enable_if<std::is_base_of<Book, typename Container::value_type::member_type>::value, double>::type
        u(const double &money, const Container &books) const {
            double u = money - penalty(books.size());
            for (auto &book : books) {
                u += uBook(book);
            }
            return u;
        }

        /** Same as the above, but takes a container of eris_id_t values of books rather than
         * SharedMember<Book> objects.
         */
        template <typename Container>
        typename std::enable_if<std::is_same<eris::eris_id_t, typename Container::value_type>::value, double>::type
        u(const double &money, const Container &books) const {
            double u = money - penalty(books.size());
            for (auto &book_id : books) {
                u += uBook(simGood<Book>(book_id));
            }
            return u;
        }

        /** Returns the utility (or potential utility) of the agent.  This is updating during the
         * intra-apply phases; before that runs, it returns the utility in the previous period.
         */
        const double& u() const;

        /** Returns the lifetime utility of the agent.  This is updated at the same time as u(),
         * adding current period utility as soon as it is finalized (during intra-apply).
         */
        const double& uLifetime() const;

        /** Returns/accesses the library of books owned by this Reader.
         *
         * \returns map where keys are the eris_id_t of the books and values are the realized
         * quality of the books. */
        const std::unordered_map<eris::eris_id_t, double>& library();

        /** Returns the set of Books that were purchased in the last period.
         *
         * This is updated only during the agent's intraApply() phase; it intra stages before
         * "apply" it gives the new books in the previous period.
         */
        const std::unordered_set<eris::SharedMember<Book>>& newBooks() const;

        /** Returns a list of eris_id_t of books this reader wrote, in chronological order from
         * earliest to latest. */
        const std::vector<eris::eris_id_t>& wrote() const;

        /** Returns the unpenalized utility of reading the given book.  The utility is found by
         * adding together the distance polynomial value and the quality value of the book.  The
         * distance polynomial is as returned by distPolynomial(), evaluated at the distance from
         * the Reader to the given Book.  The quality comes from calling the quality() method.  If
         * the calculated overall value is less than 0, 0 is returned.
         *
         * If the book is already in the reader's library, this returns the realized utility.
         * Otherwise, this method returns an estimate.  See quality() for details.
         *
         * Note that this method isn't enough to evaluate the utility of adding an additional book:
         * calling code must also ensure that the book being added is not already contained in the
         * reader's library.
         *
         * \param b the book.
         *
         * \returns the estimated utility from reading the given book, if not in the reader's
         * library; the realized utility from reading the given book if in the reader's library.
         */
        virtual double uBook(const eris::SharedMember<Book> &b) const;

        /** Returns the quality of a given book.  If the given book is already in the user's
         * library, this returns a realized quality value; otherwise it returns a predicted quality
         * value based on the reader's prior.  This quantity must be non-stochastic (calling it
         * multiple times without the library or prior having changed should return the same value).
         */
        virtual double quality(const eris::SharedMember<Book> &b) const;

        /** Utility penalty from reading `books` books.  Must be (non-strictly) increasing and
         * strictly positive, and must return 0 for 0 books.  Uses the current penalty polynomial
         * evaluated at `books`.
         */
        virtual double penalty(unsigned long books) const;

        /** Sets the coefficients for the polynomial used in uBook(double).  Coefficient \f$c_i\f$
         * is the coefficient multiplying \f$d^i\f$ in the resulting polynomial.  For example,
         * coefficients(std::vector<double>{{3, 1, -1}}) sets the coefficients for distance utility
         * of \f$-d^2 + d + 3\f$.  The polynomial is permitted to produce negative values: these will
         * be treated as 0 in uBook().
         *
         * Some rudimentary safety checks are performed to check common ways the given polynomial
         * might be invalid (note that these checks are skipped if the indicated element of `coef`
         * doesn't exist):
         * - all coef[i] values must be finite.
         * - `coef[0] > 0` must be true.  If it isn't, the decreasing requirement of uBook means
         *   that the reader gets no utility ever from books (since \f$f(x > 0) < f(0) \leq 0\f$).
         * - The first non-zero coefficient (not counting `coef[0]`) must be negative (typically
         *   this means `coef[1] <= 0`).  If it isn't, the polynomial is increasing for values of
         *   `d` close to 0 (more technically, the derivative of the polynomial is strictly positive
         *   at 0).
         * - the last non-zero coefficient (not counting `coef[0]`) must be negative.  If it isn't,
         *   the limit of the polynomial is positive infinity and thus must (eventually) be upward
         *   sloping.
         * - the polynomial must be decreasing when evaluated at 0, 0.000001, 0.001, 0.01, 0.1, 1,
         *   10, 100, 1000, 1000000.  These values are arbitrary, and of course this is no guarantee
         *   that the polynomial is decreasing between and outside these values: it only provides a
         *   safety check.
         *
         * \param coef the vector of polynomial coefficients.
         *
         * Note that the coefficients set here might not be used by a subclass that overrides
         * uBook().
         */
        void uPolynomial(std::vector<double> coef);

        /** Accesses the vector of coefficients for the uBook() polynomial.
         *
         * \sa uPolynomial(std::vector<double>)
         */
        const std::vector<double>& uPolynomial() const;

        /** Evaluates the given polynomial at the value `x`.
         *
         * \param x the value at which to evaluate the polynomial.
         * \param coefficients the coefficients at which to evaluate the polynomial.
         * `coefficients[0]` is the constant, `coefficients[i]` applies to the \f$x^i\f$ term.
         */
        static double evalPolynomial(const double &x, const std::vector<double> &polynomial);

        /** The default polynomial coefficients.  The default is the straight line \f$4-d\f$ (that
         * is, the vector {4.0, -1.0}).
         */
        static const std::vector<double> default_polynomial;

        /** Sets the coefficients for the penalty polynomial used in `penalty(n)`.
         * Coefficient \f$c_i\f$ is the coefficient multiplying \f$n^i\f$ in the resulting
         * polynomial.  For example, coefficients(std::vector<double>{{0, 1, 0.25}}) sets the
         * coefficients for a book count penalty of \f$\frac{1}{4}n^2 + n\f$.  The polynomial is
         * required to evaluate to be increasing in the non-negative integers \f$n\f$.
         *
         * Some rudimentary safety checks are performed to check common ways the given polynomial
         * might be invalid:
         * - f(0) <= f(1)
         * - f(1) <= f(2)
         * - ...
         * - f(9) <= f(10)
         * - f(10) <= f(20)
         * - f(20) <= f(30)
         * - ...
         * - f(90) <= f(100)
         * - f(100) <= f(1000)
         * - f(1000) <= f(1000000)
         *
         * \param coef the vector of polynomial coefficients.
         * \throw std::domain_error if the polynomial appears inadmissable by failure of one of the
         * above safety checks.
         */
        void penaltyPolynomial(std::vector<double> coef);

        /** Returns the coefficients for the penalty() polynomial.
         *
         * \sa penaltyPolynomial(std::vector<double>)
         */
        const std::vector<double>& penaltyPolynomial() const;

        /** The default penality polynomial coefficients.  The default is the function
         * \f$\frac{b^2}{4}\f$, which has first few values (beginning at 0 books): (0, 0.25, 1,
         * 2.25, 4)
         */
        static const std::vector<double> default_penalty_polynomial;

        /** In-between periods, the reader randomly creates a book. */
        void interApply() override;

        /** When advancing a period, the reader takes a random step. */
        void interAdvance() override;

        /** The probability that the reader will write a book in between periods. Defaults to 0.001. */
        double writer_prob{0.001};

        /** The standard deviation of a written book.  A written book will be located at the
         * reader's position plus independent draws from \f$N(0, s)\f$ in each dimension, where
         * \f$s\f$ is the value of this parameter.
         *
         * Defaults to 0.5.
         */
        double writer_book_sd{0.5};

        /// Updates the book market cache with any new books
        void intraInitialize() override;

        /** Performs intra-period optimization to pick an optimal set of books to buy this period.
         *
         * The procedure works by calculating \f$uBook(b)-p_b\f$ for all eligible books \f$b\f$, then
         * considers adding the book with the highest \f$u-p\f$ value to the current set of new
         * purchases (and subtracing p from available cash).  If u() is higher with the book than without,
         * the purchased is "locked in" and the next-highest book is considered.
         *
         * Note on eligible books: books are only considered as candidates to be added if three
         * conditions hold:
         * - the book is not already in the user's library (from past purchases)
         * - the book is not already in the set of new purchases
         * - the book's purchase price does not exceed the user's available resources
         *
         * The last point has a potential misoptimization: consider a reader with income of 2.5
         * considering three books A, B, and C with utility values of 5, 4, and 3, with prices of 2,
         * 1.5, and 1, respectively.  Assume penalties are 0.  The consumer's utility is maximized
         * by buying books B and C (u = 0 + 4-1.5 + 3-1 = 4.5), but the behaviour above will lead the
         * consumer to buy A (u = 0.5 + 5-2 = 3.5), after which neither B nor C is affordable.
         */
        void intraOptimize() override;
        /** Completes the book purchases calculated and reserved in intraOptimize(). */
        void intraApply() override;
        /** Resets the book purchases calculated in intraOptimize(). */
        void intraReset() override;

    private:
        std::vector<double> u_poly_ = Reader::default_polynomial;
        std::vector<double> pen_poly_ = Reader::default_penalty_polynomial;
        /// Map of books owned to realized quality of those books:
        std::unordered_map<eris::eris_id_t, double> library_;
        /// Reservations of books being purchased
        std::forward_list<eris::Market::Reservation> reservations_;
        /// set of books associated with reservations_
        std::unordered_set<eris::SharedMember<Book>> reserved_books_;
        /// Books purchased in the just-finished period
        std::unordered_set<eris::SharedMember<Book>> new_books_;
        /// Books written by this reader, in order from earliest to latest
        std::vector<eris::eris_id_t> wrote_;
        /// Cache of the set of book markets available
        std::unordered_set<eris::SharedMember<BookMarket>> book_cache_;

        // Quality prior coefficients.  Key is the author id, values are the author-specific
        // coefficients.  0 is for the global (non-author-specific) parameters.
        //
        // FIXME: this just returns the mean of a chi-squared-1, which is temporary for as long as a
        // chi-squared-1 is the quality distribution in interApply, but needs to be changed when
        // that changes.
        std::unordered_map<eris::eris_id_t, std::array<double, 3>> q_coef_{{0, {{1, 0, 0}}}};

        double u_curr_, u_lifetime_;
};

}
