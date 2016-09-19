#pragma once
#include <eris/Optimize.hpp>
#include <eris/WrappedPositional.hpp>
#include <eris/Agent.hpp>
#include <eris/Market.hpp>
#include "creativity/Book.hpp"
#include "creativity/BookCopy.hpp" // IWYU pragma: keep
#include "creativity/belief/ProfitStream.hpp"
#include "creativity/belief/Demand.hpp"
#include <vector>
#include <unordered_set>
#include <forward_list>
#include <functional>
#include <map>
#include <unordered_map>
#include <set>
#include <type_traits>

namespace creativity {

class Creativity;

/// Namespace for reader belief classes
namespace belief { class Profit; }

/** A Reader is an agent that both consumes previously unread books and potentially writes new books
 * and sells copies of those books.  The Reader's utility is determined by books read and an outside
 * option; in particular it is quasilinear of the form:
 *
 * \f$u(m, b_1, b_2, \ldots, b_n) = m + \sum_{i=1}^{n} f(b_i) - p(n)\f$
 *
 * where:
 *     - \f$r\f$ is the Reader's position
 *     - \f$b_1, \ldots, b_n\f$ are the \f$n\f$ books the agent buys this period
 *     - \f$f(b)\f$ is a function that maps a book into a utility value using the book's quality
 *       minus an increasing function of the distance between the reader's position and the book's
 *       position.  It returns 0 for any books that the reader already owns.
 *     - \f$p(n)\f$ is a penalty function that increases in the number of books read
 *     - \f$m\f$ is the quantity of non-book spending the reader engages in which delivers a
 *     constant marginal utility of 1.
 *
 * Behaviour
 * =========
 *
 * Reading books
 * -------------
 * A reader optimizes in each period by looking at all Books currently available for sale,
 * considering their quality, price, and distance then deciding which books to buy.  Books are only
 * purchased once, and at most one copy is purchased by any reader.  The utility gain of a book is
 * incurred immediately.
 *
 * Mean book quality is observable, but readers receive a subjective actual quality drawn from a
 * normal distribution with the observed mean and fixed standard deviation.  This subjective quality
 * is realized only when readers obtain a copy of the book.
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
 * Lifetime profit
 * ---------------
 * Readers have a belief about the lifetime profitability of a single book which depends on various
 * author attributes and the book quality.  This belief is updated based on all books this Reader
 * has obtained plus his own previous books (if any).  This is not *directly* used for prediction:
 * instead prediction uses the extrapolated beliefs, below, which incorporate profit predictions for
 * still-on-market books.
 * 
 * See creativity::belief::Profit for details.
 *
 * Lifetime profit (extrapolated)
 * ------------------------------
 * This is identical to the lifetime profit belief, above, but also includes predicted profit levels
 * of books that are still on the market (and hence have unknown final profit).  This belief is
 * discarded when advancing to the next period: the extrapolated belief is rebuilt from the next
 * period's updated profit belief.
 *
 * Profit stream beliefs
 * ---------------------
 * The reader has beliefs for books of various ages.  These stream beliefs are used to supplement
 * the lifetime profit with predictions for books still on the market in the extrapolated belief.
 *
 * See creativity::belief::ProfitStream for details.
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
 *   - If there is no in-progress book, using updated belief, the reader decides whether to create,
 *     and if so, the effort to expend in creation.
 *   - for all authored books still on the market, decides whether to remove from the market, or to
 *     keep on the market.  If the latter, also chooses the price for the upcoming period.
 * - inter-apply:
 *   - Receives external (non-book) income, less any policy taxes.
 *   - If book creation initiated in previous periods is now complete, the book is created.
 *   - If creation was decided upon in inter-optimize, pays the decided-upon effort level and fixed
 *     creation cost from income, and creates a book
 *   - If creation_time is 0, the book appears immediately
 *   - The author incurs the fixed cost for all books staying on the market (including a new one, if
 *     just created).
 *   - Takes a step of in a random direction (if reader_step_sd is positive)
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
 * - intra-finish:
 *   - (via CopyrightPolice) If the 'catch-pirates' policy is active, people are randomly caught and
 *     incur the fine (subtracted from utility).  The opportunity component of the fine disappears;
 *     the rest of the fine is redistributed to the authors of the pirated works.
 *   - (via BookMarket) The author receives all of the period's proceeds from the period.
 */
class Reader : public eris::WrappedPositional<eris::Agent>,
    public virtual eris::interopt::Begin,
    public virtual eris::interopt::OptApply,
    public virtual eris::intraopt::Initialize,
    public virtual eris::intraopt::OptApplyReset
{
    public:
        Reader() = delete; ///< Not default constructible

        /** Constructs a new Reader at the given position.
         *
         * Beliefs start off non-informative.
         *
         * \param creativity the Creativity object owning the simulation this reader is created in
         * \param pos the initial position of the reader
         */
        Reader(
                Creativity &creativity,
                const eris::Position &pos
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
            double u = money - numBooksPenalty(books.size());
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
            double u = money - numBooksPenalty(books.size());
            for (auto &book_id : books) {
                u += uBook(simGood<Book>(book_id));
            }
            return u;
        }

        /** Similar to the above, but operators on an unordered_map of book -> bookcopy-reference
         * pairs, using the quality from the BookCopy object instead of looking it up.
         */
        double u(double money, const std::unordered_map<eris::SharedMember<Book>, std::reference_wrapper<BookCopy>> &books) const;

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
         * The reader's library() contains all books purchased, pirated, and written by this reader.
         *
         * \returns map where keys are the books and values are the reader-specific BookCopy values.
         */
        const std::unordered_map<eris::SharedMember<Book>, BookCopy>& library() const;

        /** Returns the Books that were obtained in the last period, whether purchased or obtained
         * through sharing.  Books newly authored by this reader are not included.
         *
         * This is updated only during the agent's intraApply() phase; it intra stages before
         * "apply" it gives the new books in the previous period; after intraApply it contains new
         * books for the current period.
         */
        const std::unordered_map<eris::SharedMember<Book>, std::reference_wrapper<BookCopy>>& newBooks() const;

        /** Returns the set of Books that this reader wrote, sorted by id (which is also sorted by
         * creation date of the book, from oldest to newest).  These books are included in the
         * reader's library but aren't available for sharing with friends.
         */
        const std::set<eris::SharedMember<Book>>& wrote() const;

        /** Returns the unpenalized utility of reading the given book.  The utility is found by
         * subtracting the distance polynomial value from the quality value of the book.  The
         * distance polynomial is as returned by distancePenaltyPolynomial(), evaluated at the
         * distance from the Reader to the given Book.  The quality comes from calling the quality()
         * method.
         *
         * If the book is already in the reader's library, this returns the realized utility.
         * Otherwise, this method uses the mean of the quality distribution of the book.
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

        /** Just like uBook, but uses the given quality value instead of looking it up. */
        virtual double uBook(const eris::SharedMember<Book> &b, double quality) const;

        /** Returns the quality of a given book.  If the given book is already in the user's
         * library, this returns the realized quality value determined when the book was added;
         * otherwise it returns the mean of the quality distribution for the book.
         */
        virtual double quality(const eris::SharedMember<Book> &b) const;

        /** Sets the coefficients for the polynomial used in uBook(double) to calculate the penalty
         * for book distance.  Coefficient \f$c_i\f$ is the coefficient multiplying \f$d^i\f$ in the
         * resulting polynomial.  For example, `coefficients(std::vector<double>{{0, 1, 0.5}})` sets
         * the coefficients for distance penalty to of \f$d + \frac{d^2}{2}\f$.  The polynomial
         * should be increasing in \f$d\f$ for \f$d \geq 0\f$; this is checked for a few sample `d`
         * values to catch common problems, but the checking is not exhaustive.
         *
         * The first coefficient (i.e.\ the constant) is typically 0: readers receive no penalty for
         * reading a book at their most preferred location.
         *
         * Some rudimentary safety checks are performed to check common ways the given polynomial
         * might be invalid (note that these checks are skipped if the indicated element of `coef`
         * doesn't exist):
         * - all coef[i] values must be finite.
         * - the first non-zero coefficient after the constant must be non-negative (if not, the
         *   polynomial will be decreasing in the neighbourhood of \f$d=0\f$ since the first lowest
         *   order, non-negative, non-constant term dominates in values very close to 0).
         * - the last, non-zero coefficient must be positive (if not, the polynomial will
         *   (eventually) become negative.
         * - f(0) <= f(1e-100) <= f(1e-10) <= f(.1) <= f(1) <= f(10) <= f(100) <= f(1e10) must all hold (this
         *   is a rudimentary test that the function is actually increasing).
         *
         * \param coef the vector of polynomial coefficients.
         *
         * Note that the coefficients set here might not be used by a subclass that overrides
         * uBook().
         */
        void distancePenaltyPolynomial(std::vector<double> coef);

        /** Accesses the vector of coefficients for the uBook() polynomial.
         *
         * \sa uPolynomial(std::vector<double>)
         */
        const std::vector<double>& distancePenaltyPolynomial() const;

        /** The default distance penalty polynomial coefficients.  The default is {0, 0, 3}.
         */
        static constexpr std::initializer_list<double> default_distance_penalty_polynomial{0, 0, 3};

        /** Returns the distance penalty by evaluating the distancePenaltyPolynomial() at the given
         * distance value.
         */
        double distancePenalty(double distance) const;

        /** Sets the coefficients for the penalty polynomial for the number of books read used in
         * `numBooksPenalty(n)`.  Coefficient \f$c_i\f$ is the coefficient multiplying \f$n^i\f$ in
         * the resulting polynomial.  For example, coefficients(std::vector<double>{{0, 1, 0.25}})
         * sets the coefficients for a book count penalty of \f$\frac{1}{4}n^2 + n\f$.  The
         * polynomial is required to evaluate to be increasing in the non-negative integers \f$n\f$.
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
        void numBooksPenaltyPolynomial(std::vector<double> coef);

        /** Returns the coefficients for the numBooksPenalty() polynomial.
         *
         * \sa penaltyPolynomial(std::vector<double>)
         */
        const std::vector<double>& numBooksPenaltyPolynomial() const;

        /** The default penality polynomial coefficients.  The default is the function
         * \f$b^2\f$.
         */
        static constexpr std::initializer_list<double> default_num_books_penalty_polynomial{0, 0, 1};

        /** Utility penalty from reading `books` books in the same period.  The basic idea behind
         * this penalty is that reading more books has a larger opportunity cost (or equivalently,
         * books read in the same period have decreasing marginal utility).
         *
         * Must be (non-strictly) increasing and strictly positive, and must return 0 for 0 books.
         *
         * By default the current penalty polynomial is evaluated at `books`.
         */
        virtual double numBooksPenalty(unsigned long books) const;

        /** The reader's creation function shape coefficient.  This reader can exhert effort \f$\ell
         * \geq 0\f$ to create a book of quality \f$q(\ell) = \alpha \frac{(\ell+1)^\beta -
         * 1}{\beta}\f$, where \f$\beta\f$ is this value.
         *
         * Note that this function is designed to have both domain and range of the non-negative
         * rationals, and will always have \f$q(0) = 0\f$.
         *
         * Note that \f$\beta = 0\f$ is handled specially as \f$q(\ell) = \alpha ln(\ell+1)\f$ which
         * holds mathematically (by L'Hôpital's Rule). (Without this special handling, evaluating the
         * above numerically would result in a NaN value.)
         *
         * The default value is \f$\beta = 0\f$, yielding a logarithmic quality/effort relationship.
         *
         * The value of this parameter should be strictly less than 1 to maintain the concavity of
         * the function (at 1, the function becomes linear, and above 1, strictly convex).  Values
         * less than 0 are permitted as well, but introduce a horizontal asymptote that may or may
         * not be desirable: in particular they will result in functions that start at 0 and are
         * bounded by a finite asymptote at \f$\frac{\alpha}{-\beta} > 0\f$.
         */
        double creation_shape = 0.0;

        /** The reader's creation function scale coefficient, which is \f$\alpha\f$ in the equation
         * given in the documentation for `creation_shape`.
         *
         * \sa creation_shape
         *
         * The default value is \f$\alpha = 10\f$.  Changing this value across readers of a
         * simulation allows readers to differ in ability while maintaining the same functional
         * form.
         *
         * The specified value must be non-negative, for obvious reasons.  Specifying 0 is
         * technically permitted, but results in a reader always producing 0 quality works,
         * regardless of effort (and so will put no effort at all into creation).
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

        /** Returns the effort required of the reader to create a book of the given quality.  If the
         * reader cannot produce a book of the given quality at any effort level, infinity is
         * returned.
         *
         * \sa creation_shape
         * \sa creation_scale
         * \sa creationQuality
         */
        double creationEffort(double quality) const;

        /** Returns the resulting quality of a reader with the given `scale` and `shape` parameters
         * who expends effort `effort` to create a book.
         *
         * This static version of the method is called by the non-static method.
         *
         * \sa creation_shape
         * \sa creation_scale
         * \sa creationQuality
         */
        static double creationQuality(double shape, double scale, double effort);

        /** Returns the required effort of a reader with the given `scale` and `shape` parameters
         * to create a book of quality `quality`.
         *
         * This static version of the method is used by the non-static method.
         *
         * \sa creation_shape
         * \sa creation_scale
         * \sa creationEffort
         */
        static double creationEffort(double shape, double scale, double effort);

        /// Read-only access to this reader's profit belief
        const belief::Profit& profitBelief() const;
        /// Read-only access to this reader's profit belief with profit stream extrapolations
        const belief::Profit& profitExtrapBelief() const;
        /** Returns true if profitExtrapBelief is actually just profitBelief (i.e. there were no
         * profit extrapolations to incorporated this period).
         */
        bool profitExtrapBeliefDiffers() const;
        /** Read-only access to this reader's profit stream belief for books with age `age`.  The
         * returned object will be a belief::ProfitStream model with between 1 and `age` parameters.
         * When a model for the requested age is not available the model with the highest age less
         * than `age` is returned.
         *
         * If there are no models at all, a noninformative one for `age=1` is created and returned.
         *
         * \param age the requested age, which must be at least 1.
         * \param usable optional parameter that, when true (defaults to false), specifies that only
         * usable beliefs (as determined by usableBelief()) should be considered.  Note that if
         * there are no usable beliefs at all, the noninformative (or insufficiently informed)
         * `age=1` is still returned (and created, if necessary).
         *
         * Note that the returned model will have, at most, `age` parameters (that is, `model.K() <=
         * age` will always be true).
         */
        const belief::ProfitStream& profitStreamBelief(unsigned int age, bool usable = false) const;

        /** Returns true if the given belief is usable.  Specifically, to be usable it must not be a
         * noninformative() belief, and must have at least `K` plus `belief_min_n_less_k`
         * observations.  Note that this doesn't guarantee that the belief can be successfully drawn
         * from, i.e. because it has linear constraints that are too far in the tail to meaningfully
         * satisfy.
         */
        bool usableBelief(const eris::learning::BayesianLinear &model) const;

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

        /** Returns the cost (or expected cost) of obtaining a work through sharing.  The current
         * implementation of this method simply returns the `cost_piracy` parameter of the
         * Creativity object (i.e. the copying cost is the same for everyone), but this may change
         * in the future (e.g. to incorporate beliefs about the penalties of being caught).
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

        /** Called to alter the reader's utility for the current period.  This may only be called
         * after the intraApply stage has completed (at default priority), i.e. typically in a
         * late-priority intraApply or in an intraFinish.
         *
         * This is primarily intended at policy mechanisms like CopyrightPolice who need a way to
         * fine the reader; one can imagine that the utility hit is a removal of wealth, though
         * wealth isn't actually modelled in this simulation.
         */
        void alterUtility(double amount);

        /** At the beginning of the period transition, the reader takes a random step. */
        void interBegin() override;

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
        Creativity &creativity_;

        std::shared_ptr<belief::Profit> profit_belief_, ///< Belief about lifetime book profits
            profit_belief_extrap_; ///< Beliefs about lifetime book profits using profit stream expectations
        belief::Demand demand_belief_; ///< Belief about per-period demand
        /** Profit stream beliefs for books on market for various lengths of time.  E.g.
         * `profit_stream_beliefs_[3]` is the model of future profits for books that stayed on the
         * market for at least 3 periods.  Used to build profit_belief_extrap_.
         */
        std::map<unsigned int, belief::ProfitStream> profit_stream_beliefs_;

        /** The set of friends of this reader. */
        std::unordered_set<eris::SharedMember<Reader>> friends_;

        /** Updates all of the reader's beliefs, in the following order:
         * - per-period demand
         * - profit stream models
         * - lifetime profitability
         * - lifetime profitability with profit stream extrapolation for still-in-market books
         *
         * If the CreativitySettings indicates that existing beliefs should be weakened when used as
         * a prior, that is done.
         *
         * After updating beliefs, `library_on_market_` is updated to remove books that are no
         * longer on the market (since those will have now been incorporated into the beliefs).
         *
         * \sa updateDemandBelief
         * \sa updateProfitStreamBelief
         * \sa updateProfitBelief
         * \sa CreativitySettings.prior_weight
         * \sa CreativitySettings.prior_weight_piracy
         */
        void updateBeliefs();

        /** Updates the demand equation belief based on book sales observed in the previous period.
         *
         * This includes each book as a single observation with its age at time of purchase.
         *
         * \todo Consider adding books for each period they survive (and past periods).  In other
         * words, if the book is bought when \f$age=2\f$, also add data points for \f$age=0\f$ and
         * \f$age=1\f$, and keep doing so into the future until the book leaves the market.
         *
         * The current belief is used as a prior, suitably weakened if
         * CreativitySettings.prior_weight is not equal to 1.
         */
        void updateDemandBelief();

        /** Updates the lifetime profit equation belief based on observations from the previous
         * period.
         *
         * The current belief is used as a prior, suitably weakened if
         * CreativitySettings.prior_weight is not equal to 1.
         *
         * This also updates the lifetime profit extrapolation belief, using the just-updated profit
         * belief as (non-weakened) prior plus extrapolations (via the profit stream belief) for
         * books that are still on the market.  The previous extrapolated profit belief is
         * discarded.
         */
        void updateProfitBelief();

        /** Updates the profit stream equation belief based on observations from the previous
         * period.
         *
         * The current belief is used as a prior, suitably weakened if
         * CreativitySettings.prior_weight is not equal to 1.
         *
         * This method is currently disabled; attempting to call it will result in a
         * std::runtime_error exception.
         */
        void updateProfitStreamBelief();

        /** Called from Book when the Book gets added to the simulation to let the author record the
         * new book appropriately.  This can't be done at the instance the book is created because
         * it doesn't have an id yet (it won't be actually added until the end of the in-progress
         * optimization round), but needs to have an id to be put into containers.
         */
        void registerAuthoredBook(eris::SharedMember<Book> book);
        /** Called from Book::setMarket when one of this author's Books has changed market status
         * (either gaining a new private or public market, or losing its current market).
         */
        void registerMarketUpdate(eris::SharedMember<Book> book);
#ifndef DOXYGEN_SHOULD_SEE_THIS
        // The protected register* methods are called from Book:
        friend void Book::added();
        friend void Book::setMarket(eris::SharedMember<BookMarket>);
        friend void Book::weakDepRemoved(eris::SharedMember<eris::Member>, eris::eris_id_t);
#endif

    private:
        std::vector<double> dist_penalty_poly_;
        std::vector<double> nbooks_penalty_poly_;
        /// Map of books owned to properties of those books (in BookCopy container objects):
        std::unordered_map<eris::SharedMember<Book>, BookCopy> library_;
        /// Books obtained in the just-finished period
        std::unordered_map<eris::SharedMember<Book>, std::reference_wrapper<BookCopy>> library_new_;

        std::unordered_set<eris::SharedMember<Book>>
            /// Books written by this reader that are still on the market
            wrote_market_,
            /// Cache of the set of all simulation books that haven't been read yet
            book_cache_,
            /** Cache of the set of books that haven't been read and haven't been learned (because
             * they were still on the market last time profit belief was updated).  This cache is
             * cleared of off-market books in updateBeliefs().
             */
            book_cache_market_;
        /** The set of books associated with reservations_ and reserved_piracy_cost_; the status
         * indicates whether this is a market, pirated, or public book. */
        std::unordered_map<eris::SharedMember<Book>, BookCopy::Status> reserved_books_;

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

        // True in the period in which creation started, false otherwise.
        bool create_starting_ = false;
        int create_countdown_ = -1; // the number of periods remaining until the book is finished; new books can only be initiated when this is -1.
        double create_effort_ = 0, create_quality_ = 0, create_price_ = 0;
        // The author's position at the time of creation; the final position will be this plus noise
        eris::Position create_position_;

};

}
