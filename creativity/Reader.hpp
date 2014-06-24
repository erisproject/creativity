#pragma once
#include <eris/WrappedPositional.hpp>
#include <eris/agent/AssetAgent.hpp>
#include <eris/Market.hpp>
#include <string>
#include <vector>
#include <unordered_set>
#include <boost/any.hpp>
#include <Eigen/Core>
#include "belief/Profit.hpp"
#include "belief/ProfitStream.hpp"
#include "belief/Demand.hpp"
#include "belief/Quality.hpp"

namespace creativity {

class Book; // forward declaration
class BookMarket;

/** A Reader is an agent that both consumes previously unread books and potentially writes new books
 * and sells copies of those books.  The Reader's utility is determined by books read and an outside
 * option; in particular it is quasilinear of the form:
 *
 * \f$u(m, x_0, x_1, \ldots, x_{n-1}) = m + \sum_{i=1}^{n} f(b_i) - p(n)\f$
 *
 * where:
 *     - \f$r\f$ is the Reader's position
 *     - \f$b_1, \ldots, b_n\f$ are the positions of the \f$n\f$ books the agent buys this period
 *     - \f$f(b)\f$ is a function that maps a book into a utility value using a decreasing function
 *     of the distance between the reader's position and the book's position.  It returns 0 for any
 *     books that the reader already owns.
 *     - \f$p(n)\f$ is a penalty function that increases in the number of books read
 *     - \f$m\f$ is the quantity of non-book spending the reader engages in which delivers a
 *     constant marginal utility of 1.
 *
 * Behaviour
 * =========
 *
 * Reading books
 * -------------
 * A reader optimizes in each period by looking at all Books currently available for sale, assessing
 * their quality, considering their price, then deciding which books to buy.  Books are only
 * purchased once, and at most one copy is purchased by any reader.  The utility gain of a book is
 * incurred immediately.
 *
 * Book quality is not observable but book authorship is.  Readers have beliefs about book quality
 * of previously known and previously unknown authors with which they can estimate quality.  Once
 * purchased, readers obtain a random, permanent quality draw (which is based on the actual quality)
 * as their own private assessment of the book's quality.
 *
 * Writing books
 * -------------
 * "Readers" have an innate ability to write a book by using effort (measured in terms of foregone
 * income), where more effort creates a higher quality book.  This potential author has beliefs
 * (more details below) about the profitability of authorship, based on observed book sales, with
 * which he makes the decision to create or not, deciding the quality at which to create in the
 * process.
 *
 * Having created a book (or many books), the author determines a price for each of his books in
 * each period, based on a belief of the demand structure for a book.  Once the believed
 * profitability of a book becomes negative, the author permanently withdraws the book from the
 * market.
 *
 * Beliefs
 * =======
 * Reader behaviour is governed by the following beliefs which are updated over time.
 *
 * Book quality
 * ----------------------
 * When considering a book that has not been read (previously read books are not considered at all),
 * a Reader predicts quality using the model:
 *
 *     FIXME
 *     (include "known" dummy and "known" x "mean quality" terms, and perhaps "known" x "numBooks")
 *
 * to predict a book's quality.  This belief is updated whenever the reader reads any books.
 *
 * Lifetime profit
 * ---------------
 * Readers have a belief about the lifetime profitability of a single book which depends on various
 * author attributes and the book quality.  This belief is updated based on all books this Reader
 * has obtained plus his own previous books (if any).
 * 
 * See creativity::belief::Profit for details.
 *
 * Per-period demand
 * -----------------
 * Authors set a fixed price for a book in each period.  In order to do so, they have a belief about
 * the demand curve in the economy, which is updated using the price, quality, and attributes of
 * observed books.
 *
 * See creativity::belief::Demand for details.
 *
 * Eris optimization overview
 * ==========================
 * In terms of Eris stages, the decision processes are as follows:
 *
 * - inter-optimize:
 *   - updates beliefs based on previous period activity
 *   - using updated belief, decides whether to create, and if so, the effort to expend in creation
 *   - for all authored books still on the market, decides whether to remove from the market, or to
 *     keep on the market.  If the latter, also chooses the price for the upcoming period.
 * - inter-apply:
 *   - Receives external (non-book) income
 *   - If creation was decided upon in inter-optimize, gives up the decided-upon effort level in
 *     income and creates a book
 *   - The author incurs the fixed cost for all books staying on the market (including a new one, if
 *     just created).
 *   - Takes a step of distance N(0,0.25) in a random direction
 * - intra-initialize:
 *   - (internal optimization; no visible effect).
 * - intra-optimize:
 *   - Decides on a set of unowned books to purchase based on expected utility gains.
 * - intra-reset:
 *   - Resets unowned book decision variables to perform reoptimization.  Currently not actually
 *     invoked as nothing in the simulation performs intra-reoptimization.
 * - intra-apply:
 *   - Buys the books decided upon, if any.  Spends all leftover income on the numeraire good.
 * - intra-finish (via BookMarket):
 *   - The author receives all of the periods proceeds from the period.
 */
class Reader : public eris::WrappedPositional<eris::agent::AssetAgent>,
    public virtual eris::interopt::OptApply,
    public virtual eris::intraopt::Initialize,
    public virtual eris::intraopt::OptApplyReset
{
    public:
        /** Constructor takes the reader position and the (wrapping) positional boundaries, new
         * model objects for lifetime profit and per-period demand, and the fixed and per-unit costs
         * of copies.
         *
         * \param pos the initial position of the reader
         * \param b1 a vertex of the wrapping boundary box for the reader
         * \param b2 the vertex opposite `b1` of the wrapping boundary box for the reader
         * \param demand a per-period demand belief object
         * \param profit a lifetime profit belief object
         * \param quality a quality belief object
         * \param cFixed the fixed cost of keeping a book on the market
         * \param cUnit the per-unit cost of producing copies of a book
         * \param income the per-period income of the agent
         */
        Reader(
                const eris::Position &pos, const eris::Position &b1, const eris::Position &b2,
                belief::Demand &&demand, belief::Profit &&profit, belief::Quality &&quality,
                const double &cFixed, const double &cUnit, const double &income
              );

        /** Takes a money value and a container of SharedMember<Book> objects and returns the
         * utility for that set of books.  This is calculated as the money plus the sum of utilities
         * from each book (from uBook()) minus the penalty incurred by buying the given number of
         * books.
         *
         * \param money the money leftover after buying the provided set of books
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
         * \returns map where keys are the books and values are the realized quality of the books.
         */
        const std::unordered_map<eris::SharedMember<Book>, double>& library();

        /** Returns the set of Books that were purchased in the last period.
         *
         * This is updated only during the agent's intraApply() phase; it intra stages before
         * "apply" it gives the new books in the previous period.
         */
        const std::unordered_set<eris::SharedMember<Book>>& newBooks() const;

        /** Returns a list of books this reader wrote, in chronological order from earliest to
         * latest. */
        const std::vector<eris::SharedMember<Book>>& wrote() const;

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
        virtual double uBook(eris::SharedMember<Book> b) const;

        /** Returns the quality of a given book.  If the given book is already in the user's
         * library, this returns the realized quality value determined when the book was added;
         * otherwise it returns a predicted quality value based on the reader's quality prior.  This
         * quantity must be non-stochastic (i.e.  calling it multiple times with the same library
         * and prior should return the same value).
         */
        virtual double quality(eris::SharedMember<Book> b) const;

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
         * \param polynomial the coefficients at which to evaluate the polynomial.
         * `coefficients[0]` is the constant term, `coefficients[i]` applies to the \f$x^i\f$ term.
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
         * might be invalid; in particular, it must satisfy each of the following inequalities:
         *
         * \f[
         *   f(0) \leq f(1) \leq f(2) \leq \ldots \leq f(9) \leq f(10) \leq f(20) \leq f(30) \leq
         *   \ldots \leq f(90) \leq f(100) \leq f(1000) \leq f(1000000)
         * \f]
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

        /** The standard deviation of a written book.  A written book will be located at the
         * reader's position plus a draw from \f$N(0, s)\f$ in a random direction, where \f$s\f$ is
         * the value of this parameter.
         *
         * Defaults to 0.5.
         */
        double writer_book_sd{0.5};

        /** The standard deviation of the quality draw for books written by this author.  A new
         * reader realizes a quality drawn from a normal distribution with mean of the true quality
         * and this standard deviation, with negative values truncated to zero.  (Note that changing
         * negative values to 0 has the effect of increasing the mean and decreasing the variance
         * with this effect being small only when the value truncation is rare).
         */
        double writer_quality_sd{1.0};

        /** The reader's creation coefficients.  This reader can exhert effort \f$\ell \geq 0\f$ to
         * create a book of quality \f$q(\ell) = \alpha_0 + \alpha_1 \ell^{\alpha_2}\f$, where
         * \f$\alpha\f$ is this vector.
         *
         * The default value is \f$\alpha = (-2.0, 4.0, 0.25)\f$, which is appropriate for a
         * 2-dimensional world and probably a poor choice for any other dimensionality.
         */
        Eigen::Vector3d creation_coefs{-2.0, 4.0, 0.25};

        /** Returns the resulting quality of this reader exerting effort level `effort` to create a
         * book.
         *
         * \sa creation_coefs
         */
        double creationQuality(const double &effort) const;

        /** Called at the end of a period (from BookMarket::intraFinish) to transfer book revenue
         * earned during a period back to the author.  The author subtracts variable cost (if any)
         * and keeps the remaining profit.  Note that fixed costs have already been removed, and so
         * this should never result in a negative profit value.
         *
         * FIXME: make sure fixed costs are already removed!
         */
        void receiveProfits(eris::SharedMember<Book> book, const eris::Bundle &revenue);

        /// Read-only access to this reader's profit belief
        const belief::Profit& profitBelief();
        /// Read-only access to this reader's profit belief with profit stream extrapolations
        const belief::Profit& profitExtrapBelief();
        /** Read-only access to this reader's profit stream belief for books with age `age`.  The
         * returned object will be a belief::ProfitStream model with between 1 and `age` parameters.
         * When a model for the requested age is not available the model with the highest age less
         * than `age` is returned.
         *
         * If there are no models at all, a highly noninformative one for `age=1` is created and
         * returned.
         */
        const belief::ProfitStream& profitStreamBelief(unsigned long age);

        /// Read-only access to this reader's demand belief
        const belief::Demand& demandBelief();
        /// Read-only access to this reader's quality belief
        const belief::Quality& qualityBelief();

        /** In-between periods, the reader optimizes by:
         * - updates his beliefs based on characteristics of newly obtained books
         * - chooses whether or not to create a new book and, if so, the effort to expend on it
         * - for all existing books still on the market, the author decides to keep or remove a book
         *   from the market, and if keeping, the price is chosen
         */
        void interOptimize() override;

        /** Applies inter-period optimization that are visible to other agents:
         * - receives fixed, non-book income
         * - creates a book (if decided upon in interOptimize()), removing the effort cost from the
         *   just-received income.
         * - Incurs fixed costs for each book on the market (including a new one, if created)
         * - takes a random step of distance N(0,0.25) in a random direction.
         */
        void interApply() override;

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

    protected:
        /// Belief about lifetime book profits
        belief::Profit profit_belief_, ///< Belief about lifetime book profits
            profit_belief_extrap_; ///< Beliefs about lifetime book profits using profit stream expectations
        belief::Demand demand_belief_; ///< Belief about per-period demand
        belief::Quality quality_belief_; ///< Belief about book quality
        /** Profit stream beliefs for books on market for various lengths of time.  E.g.
         * `profit_stream_beliefs_[3]` is the model of future profits for books that stayed on the
         * market for at least 3 periods.
         */
        std::map<unsigned long, belief::ProfitStream> profit_stream_beliefs_;

        /// The different ages for profit_stream_beliefs
        static const std::vector<unsigned long> profit_stream_ages_;

        /** Updates all of the reader's beliefs, in the following order:
         * - book quality beliefs
         * - per-period demand
         * - profit stream models
         * - lifetime profitability
         * - lifetime profitability with profit stream extrapolation for still-in-market books
         *
         * \sa updateQualityBelief
         * \sa updateDemandBelief
         * \sa updateProfitStreamBelief
         * \sa updateProfitBelief
         */
        void updateBeliefs();

        /** Updates the quality model based on observations from the previous period.
         */
        void updateQualityBelief();

        /** Updates the demand equation belief based on book sales observed in the previous period.
         *
         * This includes each book as a single observation with its age at time of purchase.
         *
         * \todo Consider adding books for each period they survive (and past periods).  In other
         * words, if the book is bought when \f$age=2\f$, also add data points for \f$age=0\f$ and
         * \f$age=1\f$, and keep doing so into the future until the book leaves the market.
         */
        void updateDemandBelief();

        /** Updates the lifetime profit equation belief based on observations from the previous
         * period.
         *
         * This also updates the lifetime profit extrapolation belief, using the just-updated profit
         * belief as prior plus extrapolations (via the profit stream belief) for books that are
         * still on the market.  The previous extrapolated profit belief is discarded.
         */
        void updateProfitBelief();

        /** Updates the profit stream equation belief based on observations from the previous
         * period.
         */
        void updateProfitStreamBelief();

        /** Need a standard normal in various places. */
        std::normal_distribution<double> stdnormal;

    private:
        std::vector<double> u_poly_ = Reader::default_polynomial;
        std::vector<double> pen_poly_ = Reader::default_penalty_polynomial;
        /// Map of books owned to realized quality of those books:
        std::unordered_map<eris::SharedMember<Book>, double> library_;
        /// Books purchased in the just-finished period
        std::unordered_set<eris::SharedMember<Book>> library_new_;
        /// Set of books that are still on the market
        std::unordered_set<eris::SharedMember<Book>> library_market_;
        /// Reservations of books being purchased
        std::forward_list<eris::Market::Reservation> reservations_;
        /// set of books associated with reservations_
        std::unordered_set<eris::SharedMember<Book>> reserved_books_;
        /// Books written by this reader, in order from earliest to latest
        std::vector<eris::SharedMember<Book>> wrote_;
        /// Books written by this reader that are still on the market
        std::unordered_set<eris::SharedMember<Book>> wrote_market_;
        /// Cache of the set of books available
        std::unordered_set<eris::SharedMember<Book>> book_cache_;

        // Track current and cumulative utility:
        double u_curr_ = 0, u_lifetime_ = 0;

        // Costs for creating copies of a book
        double c_fixed_ = 0, c_unit_ = 0;

        // Per-period income
        double income_ = 0;

        // Book prices for the upcoming period.  If a book currently on the market isn't in here,
        // or has a negative price, it'll be removed from the market
        std::unordered_map<eris::SharedMember<Book>, double> new_prices_;

        // Whether to create, and the various attributes of that creation
        bool create_ = false;
        double create_effort_ = 0, create_quality_ = 0, create_price_ = 0;

};

}
