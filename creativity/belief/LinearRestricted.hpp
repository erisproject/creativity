#pragma once
#include <unordered_map>
#include <limits>
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <eris/noncopyable.hpp>
#include "creativity/belief/Linear.hpp"

namespace creativity { namespace belief {

/** Extension to the base Linear class that supports prior restrictions on parameters via Monte
 * Carlo integration that rejects restricted draws.
 *
 * Single variable restrictions can be specified via the lowerBounds() and upperBounds() methods
 * while more general linear restrictions can be specified by calling addRestriction().
 *
 * This implementation simply calls the base Linear::draw() repeatedly until the base draw() method
 * returns an admissible value.  If the restrictions eliminate a significant portion of the
 * probability space, this class will not work well (most draws will be thrown away).
 *
 * Like Linear, this imposes a natural conjugate (Normal-gamma) prior on the model.
 */
class LinearRestricted : public Linear {
    public:
#if !EIGEN_VERSION_AT_LEAST(3,3,0)
        /** Move constructor for Eigen versions before 3.3.  Eigen 3.2 and earlier don't have proper
         * move support, and the implicit ones break things, so we work around this by providing a
         * Move constructor that just calls the implicit copy constructor.  This, of course, means
         * that for old Eigen versions, almost nothing is saved by moving since we actually copy.
         *
         * Eigen 3.3 adds a proper move constructor, and so we don't need this: the default implicit
         * move constructor should work just fine.
         *
         * Note that LinearRestricted subclasses, so long as they aren't storing additional Eigen
         * types, can rely on their default move constructors.
         */
        LinearRestricted(LinearRestricted &&move) : LinearRestricted(move) {}
        /// Default constructor
        LinearRestricted() = default;
        /// Default copy constructor
        LinearRestricted(const LinearRestricted &copy) = default;
        /// Default copy assignment operator
        LinearRestricted& operator=(const LinearRestricted &copy) = default;
#endif

        /// Constructor inherited from Linear
        using Linear::Linear;

    protected:
        // forward declaration
        class RestrictionProxyList;

    public:

        /** Returns a subscriptable object that allows assigning coefficient upper-bounds.  Indexed
         * element `i` is the upper bound for `beta[i]`, for `i` values between 0 and K-1.
         *
         * For example, to add the restriction \f$\beta_2 \leq 5\f$:
         * 
         *     model.upperBounds()[2] = 5.0;
         *
         * NaN and infinity (positive or negative) are treated as non-binding values.  The default
         * for all parameters is std::numeric_limits<double>::quiet_NaN().  To remove a previously
         * set restriction, set its value to NaN:
         *
         *     model.upperBounds()[2] = std::numeric_limits<double>::quiet_NaN();
         *
         * Warning: no checking is performed for the viability of specified bounds: specifying
         * conflicting upper and lower bounds (e.g. \f$\beta_2 \geq 3\f$ and \f$\beta_2 \leq 2\f$)
         * will result in draw() never returning, as will nonsensical restrictions such as \f$s^2
         * \leq 0\f$.
         */
        RestrictionProxyList upperBounds();

        /** Returns a random access iterator (which can be dereferenced via subscripting) to the
         * beginning of the vector of lower bounds for the \f$\beta\f$ parameters.  Element `i` is
         * the lower bound for `beta[i]`, for `i` values between 0 and K-1, and is the upper bound
         * for `s2` for `i=K`.
         *
         * For example, to add the restriction \f$\beta_2 \geq 0\f$:
         * 
         *     model.lowerBounds()[2] = 0.0;
         *
         * NaN and infinity (positive or negative) are treated as non-binding values.  The default
         * for all parameters is std::numeric_limits<double>::quiet_NaN().  To remove a previously
         * set restriction, set its value to NaN:
         *
         *     model.lowerBounds()[2] = std::numeric_limits<double>::quiet_NaN();
         *
         * Warning: no checking is performed for the viability of specified bounds: specifying
         * conflicting upper and lower bounds (e.g. \f$\beta_2 \geq 3\f$ and \f$\beta_2 \leq 2\f$)
         * will result in draw() never returning.
         */
        RestrictionProxyList lowerBounds();

        /** Adds a restriction of the form \f$R\beta \leq r\f$, where \f$R\f$ is a 1-by-`K()` row
         * vector selecting \f$\beta\f$ elements.  For example, to add a \f$\beta_2 \in [-1, 3.5]\f$
         * restriction, either of the following two approaches can be used:
         *
         *     model.upperBounds()[2] = 3.5;
         *     model.lowerBounds()[2] = -1;
         *
         *     Eigen::RowVectorXd R = Eigen::RowVectorXd::Zero(K());
         *     R[2] = 1;
         *     model.addRestriction(R, 3.5);
         *     R[2] = -1;
         *     model.addRestriction(R, 1);
         *
         * NB: to add a greater-than-or-equal restriction, convert the restriction to a
         * less-than-or-equal one (by negating the `R` and `r` values, as in the above example).
         *
         * This method allows arbitrary linear relationships between the variables.  For example,
         * the following adds the restriction \f$\beta_1 - 3\beta_3 \geq 2.5\f$:
         *
         *     Eigen::RowVectorXd R = Eigen::RowVectorXd::Zero(K());
         *     R[1] = -1; R[3] = 3;
         *     model.addRestriction(R, -2.5);
         *
         * Example adding the restriction \f$\beta_1 \leq 10\beta_5\f$ to a model with 7
         * parameters, using Eigen's comma initialization syntax:
         *
         *     Eigen::RowVectorXd R(7);
         *     R << 0, 1, 0, 0, 0, -10, 0;
         *     model.addRestriction(R, 0);
         *
         * \throws std::logic_error if R has a length other than `K()`.
         */
        void addRestriction(const Eigen::Ref<const Eigen::RowVectorXd> &R, double r);

        /** Adds a `RB >= r` restriction.  This is simply a shortcut for calling
         *
         *     addRestriction(-R, -r)
         *
         * \sa addRestriction
         */
        void addRestrictionGE(const Eigen::Ref<const Eigen::RowVectorXd> &R, double r);

        /** Adds a set of linear restrictions of the form \f$R\beta \leq r\f$, where \f$R\f$ is a
         * l-by-`K()` matrix of parameter selection coefficients and \f$r\f$ is a l-by-1 vector of
         * value restrictions corresponding to the `l` restrictions specified in `R`.  The
         * restrictions must be satisfied for all `l` rows.
         *
         * This is equivalent to using addRestriction() l times, on each row of `R` and `r`.
         *
         * `R` and `r` must have the same number of rows.
         */
        void addRestrictions(const Eigen::Ref<const Eigen::MatrixXd> &R, const Eigen::Ref<const Eigen::VectorXd> &r);

        /** Adds a set of linear restrictions of the form \f$R\beta \geq r\f$.  This is simply a
         * shortcut for calling
         *
         *     addRestrictions(-R, -r);
         *
         * \sa addRestrictions
         */
        void addRestrictionsGE(const Eigen::Ref<const Eigen::MatrixXd> &R, const Eigen::Ref<const Eigen::VectorXd> &r);

        /** Clears all current model restrictions. */
        void clearRestrictions();

        /** Overridden to reset the cached value of mean beta draws.  Since draws produced by this
         * class are independent, this does not actually perform any draws.
         *
         * Calling discard(0) will perform no actual burn-in (even for subclasses that might
         * override to provide a burn-in), but will ensure that the next predict() call is based on
         * a fresh set of beta draws (instead of reusing the results of previous predict() draws).
         */
        virtual void discard(unsigned int burn) override;

        /** Draws a set of parameter values.  initialize() must have been called before calling this
         * method.  Returns a const reference to the `k+1` length vector where elements `0` through
         * `k-1` are the \f$\beta\f$ values and element `k` is the \f$s^2\f$ value.
         *
         * \throws std::runtime_error if, due to the imposed restrictions, at least
         * `draw_max_discards` sequential draws are discarded without finding an admissible draw.
         */
        virtual const Eigen::VectorXd& draw() override;

        /** Predicts using this model using the beta values from at least 1000 draws.
         *
         * The mean of the drawn beta values are cached and will be reused by the next call to
         * draw() unless discard() is called between calls to draw.
         *
         * This method simply calls `predict(Xi, 1000)`.
         */
        virtual double predict(const Eigen::Ref<const Eigen::RowVectorXd> &Xi) override;

        /** Predicts using this model using at least the given number of beta values.
         *
         * If the currently cached beta means contains fewer than the requested number of draws, it
         * is updated with new draws before being used.  If a previous predict() call specified a
         * larger number of draws, the beta means from the larger number of draws is used.
         *
         * As a consequence, `double y1 = predict(X, 1000); double y2 = predict(X, 1000)` will
         * return the same predictions, but `double y1 = predict(X, 1000); predict(X, 2000); double
         * y2 = predict(X, 1000)` could produce different values of `y1` and `y2`.
         *
         * \throws std::logic_error if attempting to call predict() on an empty or noninformative
         * model.  (Because of the large s2 and small n values, draws would be highly random and in
         * no way useful for prediction).
         */
        virtual double predict(const Eigen::Ref<const Eigen::RowVectorXd> &Xi, long min_draws);

        unsigned long draw_discards = 0, ///< Tracks the number of draws discarded by the most recent call to draw()
                      draw_discards_max = 1000, ///< The maximum number of discards in a single draw before raising an exception
                      draw_discards_cumulative = 0, ///< Tracks cumulative draw discards over the life of this object
                      draw_success_cumulative = 0; ///< Tracks the number of successful draws over the life of this object

        /// Accesses the restriction coefficient selection matrix (the \f$R\f$ in \f$R\beta <= r\f$).
        const Eigen::Ref<const Eigen::MatrixXd> R() const;

        /// Accesses the restriction value vector (the \f$r\f$ in \f$R\beta <= r\f$).
        const Eigen::Ref<const Eigen::VectorXd> r() const;

    protected:
        /// Creates a LinearRestricted from a Linear rvalue
        LinearRestricted(Linear &&move) : Linear(std::move(move)) {}

        /** Overridden to also reset mean_beta_draws_ to 0 (to force predict() value redrawing). */
        virtual void discardForce(unsigned int burn) override;

        /** Proxy object that converts assignments into restriction rows on the associated
         * LinearRestricted model.
         */
        class RestrictionProxy final {
            public:
                /** Adds a restriction on the referenced parameter.
                 *
                 * \sa RestrictionProxyList
                 */
                RestrictionProxy& operator=(double r);
                /** Returns true if there is a restriction of the appropriate type (upper- or
                 * lower-bound) on the given parameter.  This method works whether the restriction
                 * was added via a RestrictionProxy object, or as a single-parameter restriction
                 * given to addRestriction() (or variants).
                 *
                 * Specifically this looks for any restriction with a negative (for lower bounds) or
                 * positive (for upper bounds) selector for the requested coefficient and 0 for all
                 * other coefficient selectors.
                 */
                bool restricted() const;
                /** Returns the appropriate value restriction, if there is a single-parameter
                 * restriction of the given type.  If there is no such restriction, returns a
                 * quiet_NaN.  If there are multiple restrictions, the most binding restriction is
                 * returned (that is, 4 is returned if both `>= 3` and `>= 4` lower-bound
                 * restrictions are found).
                 *
                 * Note that restrictions added with non-unitary coefficients are handled properly;
                 * that is, if a `<=` restriction on a 4-parameter model of R=(0, -2.5, 0, 0), r=1
                 * is added, the returned lower bound restriction will be -0.4.
                 */
                operator double() const;
            private:
                /** Creates a new RestrictionProxy object that, when assigned to, adds an upper or
                 * lower-bound restriction on beta[`k`].
                 *
                 * \param the LinearRestricted object where the restriction will be added;
                 * \param k the parameter index
                 * \param upper true if this is an upper bound, false if a lower bound
                 */
                RestrictionProxy(LinearRestricted &lr, size_t k, bool upper);
                friend class RestrictionProxyList;
                LinearRestricted &lr_;
                const size_t k_;
                const bool upper_;
        };
        /** Object returned by upperBounds() and lowerBounds() that provides subscript access to
         * add individual parameter restrictions (via a RestrictionProxy object).
         */
        class RestrictionProxyList final {
            public:
                /** Creates and returns a RestrictionProxy object which can be assigned a double
                 * value to add a restriction on the given parameter to the LinearRestricted model.
                 */
                RestrictionProxy operator[](size_t k);
            private:
                RestrictionProxyList(LinearRestricted &lr, bool upper);
                friend class LinearRestricted;
                LinearRestricted &lr_;
                size_t k_;
                bool upper_;
        };

        /** Stores the coefficient selection matrix for arbitrary linear restrictions passed to
         * addRestriction() or addRestrictions(). Note that the only the first
         * `restrict_linear_size_` rows of the matrix will be set, but other rows might exist with
         * uninitialized values. */
        Eigen::MatrixXd restrict_select_;
        /** Stores the value restrictions for arbitrary linear restrictions passed to
         * addRestriction() or addRestrictions().  Note that the only the first
         * `restrict_linear_size_` values of the vector will be set, but other values might exist
         * with uninitialized values. */
        Eigen::VectorXd restrict_values_;
        /** Stores the number of arbitrary linear restrictions currently stored in
         * restrict_linear_select_ and restrict_linear_values_.
         */
        size_t restrict_size_ = 0;

        /** Called to ensure the above are set and have (at least) the required number of rows free
         * (beginning at row `restrict_linear_size_`). */
        void allocateRestrictions(size_t more);

        /// The cache of drawn beta vectors used for prediction.
        Eigen::VectorXd mean_beta_;

        /// The number of beta draws used to calculate mean_beta_
        long mean_beta_draws_ = 0;
};

}}
