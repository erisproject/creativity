#pragma once
#include <eris/Positional.hpp>
#include <eris/agent/AssetAgent.hpp>
#include <string>
#include <vector>

namespace creativity {

/** A Reader is an agent with a position whose utility is determined by books and an outside option.
 * In particular, his utility is quasilinear of the form:
 *
 * \f$u(m, x_0, x_1, \hdots, x_{n-1}) = m + \sum_{i=0}^{n-1} w_i f(\| p - x_i\|)\f$
 *
 * where \f$p\f$ is the Reader's position, $x_0, \hdots, x_{n-1}$ are the position of $n$ books the
 * agent reads this period.  \f$f(d)\f$ is assumed to be decreasing in \f$d\f$ (i.e. books that
 * are further away deliver lower utility).  The weights, \f$w_0, \hdots, w_{n-1}\f$ are also
 * decreasing (i.e. \f$w_0 \geq w_1 \geq \hdots \geq w_{n-1}\f$), reflecting that additional books
 * read in a period decrease utility more than the first book.  \f$m\f$ is the quantity of non-book
 * spending the reader engages in which delivers constant marginal utility of 1.
 *
 * The reader optimizes by looking at all Books currently available for sale, assessing their
 * quality, considering their price, then deciding which books to buy (and, implicitly, which order
 * to read them).
 *
 * The default weights (which can be accessed and adjusted by calling weights() and
 * weights(std::vector<double>)) are stored in the default_weights initializer list.
 */
class Reader : public eris::Positional<eris::agent::AssetAgent>,
    public virtual eris::interopt::Apply,
    public virtual eris::intraopt::OptApplyReset
{
    public:
        /// Inherit positional constructor
        using eris::Positional<eris::agent::AssetAgent>::Positional;

        /** Takes any type of container of numerical values of distances and constructs a utility
         * value by adding together the u_part() values for each of them.
         *
         * \param distances a standard, in-order container (for example, std::vector<double>) of numeric distances
         */
        template <typename Container, typename = typename std::enable_if<std::is_arithmetic<typename Container::value_type>::value>::type>
        double u(const Container distances) const {
            size_t i = 0;
            double ret = 0.0;
            for (auto &d : distances) ret += u_nth_book(d, i++);
            return ret;
        }

        /** Utility penalty from reading `books` books.  Must be (non-strictly) increasing and
         * strictly positive, and must return 0 for 0 books.  The default is \f$\frac{b^2}{4}\f$
         */
        virtual double penalty(unsigned int books) const;

        /** Returns the utility of a book at distance `d` which is the reader's `n+1`th book (i.e. 0
         * is the first book).
         *
         * This implementation (subclasses could modify) calculates \f$weights_n u_d(d)\f$ if `n` is less
         * than the length of `weights`, and 0 otherwise.
         */
        virtual double uNthBook(double d, unsigned long n) const;

        /** Returns the unweighted utility of a book at distance `d`.  The return value of this
         * function must be finite and non-negative for all non-negative values of `d` (including
         * 0), and moreover must be (non-strictly) decreasing in `d`.
         *
         * The default implementation returns a polynomial of `d` with coefficients as given by
         * polynomial(), or 0, if the polynomial at `d` evaluates to a negative number.
         *
         * \param d the distance to the book.  Must be non-negative.
         *
         * \throws std::domain_error if `d` is negative.
         */
        virtual double uBook(double d) const;

        /** Sets the coefficients for the polynomial coefficients used in u_book(double).
         * Coefficient \f$c_i\f$ is the coefficient multiplying \f$d^i\f$ in the resulting
         * polynomial.  For example, coefficients(std::vector<double>{{3, 1, -1}}) sets the
         * coefficients for distance utility of \f$-d^2 + d + 3.  The polynomial is permitted to
         * produce negative values: these will be treated as 0 in u_book().
         *
         * Some rudimentary safety checks are performed to check common ways the given polynomial
         * might be invalid (note that these checks are skipped if the indicated element of `coef`
         * doesn't exist):
         * - all coef[i] values must be finite.
         * - `coef[0] > 0` must be true.  If it isn't, the decreasing requirement of u_book means
         *   that the reader gets no utility ever from books (since \f$f(x > 0) < f(0) \leq 0\f$).
         * - The first non-zero coefficient (not counting `coef[0]`) must be negative (typically
         *   this means `coef[1] < 0`).  If it isn't, the polynomial is increasing for values of `d`
         *   close to 0 (more technically, the derivative of the polynomial is strictly positive at
         *   values close to 0).
         * - the last non-zero coefficient (not counting `coef[0]`) must be negative.  If it isn't,
         *   the limit of the polynomial is positive infinity and thus must (eventually) be upward
         *   sloping.
         * - the polynomial must be decreasing when evaluated at 0, 0.000001, 0.001, 0.01, 0.1, 1,
         *   10, 100, 1000, 1000000.  These values are arbitrary, and if course this is no guarantee
         *   that the polynomial is decreasing between and outside these values: it only provides a
         *   safety check.
         *
         * \param coef the vector of polynomial coefficients.
         */
        void polynomial(std::vector<double> coef);

        /** Returns the coefficients for the u_book() polynomial.
         *
         * \sa polynomial(std::vector<double>)
         */
        const std::vector<double>& polynomial() const;

        /** Evaluates the polynomial defined by polynomial() at the value `x`.
         *
         * \param x the value at which to evaluate the polynomial.
         */
        double evalPolynomial(const double &x) const;

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

        /** Sets the weights for this reader.  Weights must be (non-strictly) decreasing and
         * non-negative.  Weights for book indexes past the end of the weights vector will be
         * implicitly considered to be 0.
         *
         * \param weights an std::vector<double> of new weights to assign
         *
         * \throws std::domain_error if weights is non-decreasing or has negative values
         */
        void weights(std::vector<double> weights);

        /** Accesses the current weights for this reader. */
        const std::vector<double>& weights() const;

        /** The default weights for new Reader objects: the vector {1.0, 0.75, 0.5, 0.25} */
        static const std::vector<double> default_weights;

        /** In-between periods, the reader randomly creates a book. */
        void interApply() override;

        /** The probability that the reader will write a book in between periods. Defaults to 0.001. */
        double writer_prob{0.001};

        /** The standard deviation of a written book.  A written book will be located at the
         * reader's position plus independent draws from \f$N(0, s)\f$ in each dimension, where
         * \f$s\f$ is the value of this parameter.
         *
         * Defaults to 0.5.
         */
        double writer_book_sd{0.5};

        // Optimize: uses utility function 
        void intraOptimize() override;
        void intraApply() override;
        void intraReset() override;

    private:
        std::vector<double> polynomial_ = Reader::default_polynomial;
        std::vector<double> weights_ = Reader::default_weights;
};

}
