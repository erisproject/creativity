#pragma once
#include <eris/WrappedPositional.hpp>
#include <eris/agent/AssetAgent.hpp>
#include <eris/Market.hpp>
#include <Eigen/Core>
#include "belief/Profit.hpp"
#include "belief/ProfitStream.hpp"
#include "belief/Demand.hpp"
#include "belief/Quality.hpp"
#include <string>
#include <vector>
#include <unordered_set>

namespace creativity {

class Book;
class BookMarket;
class Creativity;

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
 *   - Decides on a set of unowned books to purchase based on expected utility gains.  If piracy has
 *     been invented, the user also considers works that can be obtained as shared copies from friends.
 * - intra-reset:
 *   - Resets unowned book decision variables to perform reoptimization.  Currently not actually
 *     invoked as nothing in the default creativity simulation triggers intra-reoptimization.
 * - intra-apply:
 *   - Buys and/or pirates the books decided upon.  Spends all leftover income on the numeraire
 *     good.
 * - intra-finish (via BookMarket):
 *   - The author receives all of the periods proceeds from the period.
 */
class Reader : public eris::WrappedPositional<eris::agent::AssetAgent>,
    public virtual eris::interopt::OptApply,
    public virtual eris::intraopt::Initialize,
    public virtual eris::intraopt::OptApplyReset
{
    public:
        Reader() = delete; ///< Not default constructible

        /** Constructor takes the reader position and various reader properties.
         *
         * ProfitStream beliefs start off with a highly non-informative prior for age=1.
         *
         * \param creativity the Creativity object owning the simulation this reader is created in
         * \param pos the initial position of the reader
         * \param demand a per-period demand belief object
         * \param profit a lifetime profit belief object
         * \param quality a quality belief object
         * \param cFixed the fixed cost of keeping a book on the market
         * \param cUnit the per-unit cost of producing copies of a book
         * \param income the per-period income of the agent
         */
        Reader(
                std::shared_ptr<Creativity> creativity,
                const eris::Position &pos,
                belief::Demand &&demand, belief::Profit &&profit, belief::Quality &&quality,
                double cFixed, double cUnit, double income
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
        u(double money, const Container &books) const {
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
        u(double money, const Container &books) const {
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
         * The reader's library() contains the books returned by libraryPurchased(),
         * libraryPirated(), and wrote().
         *
         * \returns map where keys are the books and values are the realized quality of the books.
         */
        const std::unordered_map<eris::SharedMember<Book>, double>& library() const;

        /** Returns/accesses the set of books in this reader's library that were purchased by the
         * Reader.
         *
         *
         * The reader's library() contains the books returned by libraryPurchased(),
         * libraryPirated(), and wrote().
         */
        const std::unordered_set<eris::SharedMember<Book>>& libraryPurchased() const;

        /** Returns/accesses the set of books in this reader's library that were pirated by the
         * Reader.
         *
         *
         * The reader's library() contains the books returned by libraryPurchased(),
         * libraryPirated(), and wrote().
         */
        const std::unordered_set<eris::SharedMember<Book>>& libraryPirated() const;

        /** Returns the set of Books that were obtained in the last period, whether purchased or
         * obtained through sharing.  Books newly authored by this reader are not included.
         *
         * This is updated only during the agent's intraApply() phase; it intra stages before
         * "apply" it gives the new books in the previous period.
         */
        const std::unordered_set<eris::SharedMember<Book>>& newBooks() const;

        /** Returns the set of Books that were purchased in the last period, excluding pirated
         * works.
         *
         * This is updated only during the agent's intraApply() phase; it intra stages before
         * "apply" it gives the new books in the previous period.
         */
        const std::unordered_set<eris::SharedMember<Book>>& newPurchased() const;

        /** Returns the set of Books that were obtained through sharing (not purchasing) in the last
         * period.
         *
         * This is updated only during the agent's intraApply() phase; it intra stages before
         * "apply" it gives the new books in the previous period.
         */
        const std::unordered_set<eris::SharedMember<Book>>& newPirated() const;

        /** Returns the set of Books that this reader wrote, sorted by creation date of the book.
         * These books are included in the reader's library but aren't available for sharing with
         * friends.
         *
         * The reader's library() contains the books returned by libraryPurchased(),
         * libraryPirated(), and wrote().
         */
        const std::set<eris::SharedMember<Book>>& wrote() const;

        /** Returns the unpenalized utility of reading the given book.  The utility is found by
         * adding together the distance polynomial value and the quality value of the book.  The
         * distance polynomial is as returned by distPolynomial(), evaluated at the distance from
         * the Reader to the given Book.  The quality comes from calling the quality() method.  If
         * the calculated overall value is less than 0, 0 is returned.
         *
         * If the book is already in the reader's library, this returns the realized utility.
         * Otherwise, this method returns an estimate by calling quality().
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
         * quantity may be stochastic for the initial call, but subsequent calls will return the same
         * predicted value until the quality belief is updated.
         */
        virtual double quality(eris::SharedMember<Book> b) const;

        /** Utility penalty from reading `books` books in the same period.  The basic idea behind
         * this penalty is that reading more books has a larger opportunity cost (or equivalently,
         * books read in the same period have decreasing marginal utility).
         *
         * Must be (non-strictly) increasing and strictly positive, and must return 0 for 0 books.
         *
         * By default the current penalty polynomial is evaluated at `books`.
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
        static double evalPolynomial(double x, const std::vector<double> &polynomial);

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

        /// The fixed cost of keeping a book on the market for a period.
        double cost_fixed{0};

        /** The unit cost of producing a copy of a book on the market.  Copies are produced
         * instantly as required; no preallocated number of copies is performed.
         */
        double cost_unit{0};

        /// The per-unit income the reader receives.
        double income{0};

        /** The reader's creation function shape coefficient.  This reader can exhert effort \f$\ell
         * \geq 0\f$ to create a book of quality \f$q(\ell) = \beta \frac{(\ell+1)^{1-\alpha} -
         * 1}{1 - \alpha}\f$, where \f$\alpha\f$ is this value.
         *
         * Note that \f$\alpha = 1\f$ is handled specially as \f$q(\ell) = \beta ln(\ell+1)\f$ which
         * holds mathematically (by L'Hôpital's Rule); without this special handling, evaluating the
         * above numerically would result in a NaN value.
         *
         * The default value is \f$\alpha = 1\f$, yield a logarithmic quality/effort relationship.
         *
         * The value of this parameter should be strictly greater than 0 to maintain the concavity
         * of the function.
         */
        double creation_shape = 1.0;

        /** The reader's creation function scale coefficient.  This reader can exhert effort \f$\ell
         * \geq 0\f$ to create a book of quality \f$q(\ell) = \beta \frac{(\ell+1)^{1-\alpha} -
         * 1}{1 - \alpha}\f$, where \f$\beta\f$ is this value.
         *
         * The default value is \f$beta = 10\f$.  Changing this value across readers of a simulation
         * allows readers to differ in ability while maintaining the same functional form.
         *
         * The specified value must be non-negative, for obvious reasons.  Specifying 0 works as
         * expected (the reader always produces quality 0 works, regardless of effort).
         */
        double creation_scale = 10.0;

        /** Returns the resulting quality of this reader exerting effort level `effort` to create a
         * book.
         *
         * \sa creation_shape
         * \sa creation_scale
         * \sa creationEffort
         */
        double creationQuality(double effort) const;

        /** Returns the effort required of the reader to create a book of the given quality.
         *
         * \sa creation_shape
         * \sa creation_scale
         * \sa creationQuality
         */
        double creationEffort(double quality) const;

        /** Called at the end of a period (from BookMarket::intraFinish) to transfer book revenue
         * earned during a period back to the author.  The author subtracts variable cost (if any)
         * and keeps the remaining profit.  Note that fixed costs have already been removed, and so
         * this should never result in a negative profit value.
         *
         * FIXME: make sure fixed costs are already removed!
         */
        void receiveProfits(eris::SharedMember<Book> book, const eris::Bundle &revenue);

        /// Read-only access to this reader's profit belief
        const belief::Profit& profitBelief() const;
        /// Read-only access to this reader's profit belief with profit stream extrapolations
        const belief::Profit& profitExtrapBelief() const;
        /** Read-only access to this reader's profit stream belief for books with age `age`.  The
         * returned object will be a belief::ProfitStream model with between 1 and `age` parameters.
         * When a model for the requested age is not available the model with the highest age less
         * than `age` is returned.
         *
         * If there are no models at all, a highly noninformative one for `age=1` is created and
         * returned.
         *
         * Note that the returned model will have, at most, `age` parameters (that is, `model.K() <=
         * age` will always be true).
         */
        const belief::ProfitStream& profitStreamBelief(unsigned int age) const;

        /** Returns the map of all current profit stream beliefs.  The keys of the map are the
         * minimum age for which the belief applies and the value is the actual belief with `age`
         * (i.e. the key) parameters.
         *
         * Note that the model associated with key `x` will have exactly `x` parameters.
         */
        const std::map<unsigned int, belief::ProfitStream>& profitStreamBeliefs() const;

        /** The different ages that will be considered for profit_stream_beliefs.  Note that these
         * aren't won't actually have associated beliefs until a book of the given market age is
         * observed: these are only the possible age model sizes that will be used.
         */
        static const std::vector<unsigned int> profit_stream_ages;

        /// Read-only access to this reader's demand belief
        const belief::Demand& demandBelief() const;
        /// Read-only access to this reader's quality belief
        const belief::Quality& qualityBelief() const;

        /** Returns the cost (or expected cost) of obtaining a work through sharing.  The current
         * implementation of this method simply returns the `cost_unit` parameter of the Creativity
         * object (i.e. the copying cost is the same for individuals as for the creator), but this
         * may change in the future (e.g. to incorporate beliefs about being caught).
         */
        double piracyCost() const;

        /** Read-only access to the set of friends of this reader. */
        const std::unordered_set<eris::SharedMember<Reader>>& friends() const;

        /** Adds a friend to this reader, and adds this reader as a friend of the other reader.
         *
         * If the friendship is already established, this method does nothing.
         *
         * \param new_pal the new friend
         * \param recurse true (the default) if the reciprocal friendship should be added.  Internal
         * use only; external callers should not specify this value.
         *
         * \returns true if the call established a new friendship, false if the friendship already
         * existed.
         */
        bool addFriend(eris::SharedMember<Reader> new_pal, bool recurse = true);

        /** Removes a friend from this reader, and removes this reader from the given reader's set
         * of friends.
         *
         * \param old_pal the friend to remove
         * \param recurse true (the default) if the reciprocal friendship should also be removed.
         * Internal use only; external callers should not specify this value.
         *
         * \returns true if the friendship was removed, false if the friendship did not exist.
         */
        bool removeFriend(const eris::SharedMember<Reader> &old_pal, bool recurse = true);

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
        /// The Creativity object that owns the simulation this reader belongs to
        std::shared_ptr<Creativity> creativity_;
        belief::Profit profit_belief_, ///< Belief about lifetime book profits
            profit_belief_extrap_; ///< Beliefs about lifetime book profits using profit stream expectations
        belief::Demand demand_belief_; ///< Belief about per-period demand
        belief::Quality quality_belief_; ///< Belief about book quality
        /** Profit stream beliefs for books on market for various lengths of time.  E.g.
         * `profit_stream_beliefs_[3]` is the model of future profits for books that stayed on the
         * market for at least 3 periods.
         */
        std::map<unsigned int, belief::ProfitStream> profit_stream_beliefs_;

        /** The set of friends of this reader. */
        std::unordered_set<eris::SharedMember<Reader>> friends_;

        /** Updates all of the reader's beliefs, in the following order:
         * - book quality beliefs
         * - per-period demand
         * - profit stream models
         * - lifetime profitability
         * - lifetime profitability with profit stream extrapolation for still-in-market books
         *
         * After updating beliefs, `library_unlearned_` is updated to remove books that are no
         * longer on the market (since those will have now been incorporated into the beliefs).
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
        /** Map of books to quality predictions (which are cached until the library/quality belief
         * changes). */
        mutable std::unordered_map<eris::SharedMember<Book>, double> quality_predictions_;
        std::unordered_set<eris::SharedMember<Book>>
            /// Books in `library_` that were purchased by this reader
            library_purchased_,
            /// Books in `library_` that were pirated by this reader
            library_pirated_,
            /// Books obtained in the just-finished period
            library_new_,
            /// Books *purchased* in the just-finished period
            library_new_purchased_,
            /// Books obtained through sharing in the just-finished period
            library_new_pirated_,
            /** Set of books that haven't been used for learning yet (either because they are still
             * on the market, or because they were just pirated). */
            library_unlearned_,
            /// Books written by this reader that are still on the market
            wrote_market_,
            /// Cache of the set of books that haven't been read yet
            book_cache_;
        /** The set of books associated with reservations_ and reserved_piracy_cost_; the bool is
         * true if this is a pirated book, false for a purchased book. */
        std::unordered_map<eris::SharedMember<Book>, bool> reserved_books_;

        /** Books in `library_` that were authored by this reader.  The set is sorted by book ID;
         * since IDs are monotonic, this also means books are sorted in creation order.
         */
        std::set<eris::SharedMember<Book>> wrote_;

        /// Reservations of books being purchased
        std::forward_list<eris::Market::Reservation> reservations_;
        /// Total cost of piracy for all books being pirated this period
        double reserved_piracy_cost_ = 0.0;

        // Track current and cumulative utility:
        double u_curr_ = 0, u_lifetime_ = 0;

        // Book prices for the upcoming period.  If a book currently on the market isn't in here,
        // or has a negative price, it'll be removed from the market
        std::unordered_map<eris::SharedMember<Book>, double> new_prices_;

        // Whether to create, and the various attributes of that creation
        bool create_ = false;
        double create_effort_ = 0, create_quality_ = 0, create_price_ = 0;

};

}
