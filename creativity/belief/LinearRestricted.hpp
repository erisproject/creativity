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
 * Single variable restrictions can be specified via the lowerBound(), upperBound(), and restrict()
 * methods while more general linear restrictions can be specified by calling addRestriction().
 *
 * This implementation can draw by either performing Gibbs sampling or by using rejection sampling
 * (calling Linear::draw() repeatedly until the base draw() method returns an admissible value).  By
 * default, it attempts rejection sampling first and switches to Gibbs sampling if rejection
 * sampling fails.
 *
 * Like Linear, this imposes a natural conjugate (Normal-gamma) prior on the model.
 *
 * Warning: no checking is performed for the viability of restrictions when they are added:
 * specifying impossible-to-satisfy restrictions (for example, both \f$\beta_1 + 2 \beta_2 \geq 4\f$
 * and \f$2\beta_1 + 4\beta_2 \leq 2\f$) upper and lower bounds (e.g. \f$\beta_2 \geq 3\f$ and
 * \f$\beta_2 \leq 2\f$) will result in drawRejection() always failing, and an error from
 * drawGibbs().
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
        class RestrictionProxy;
        class RestrictionIneqProxy;

    public:

        /** Returns a proxy object that allows assigning coefficient upper bounds.  Element `i` is
         * the upper bound for `beta[i]`, for `i` values between 0 and K-1.
         *
         * For example, to add the restriction \f$\beta_2 \leq 5\f$:
         * 
         *     model.upperBound(2) = 5.0;
         *
         * Warning: no checking is performed for the viability of specified bounds: specifying
         * conflicting upper and lower bounds (e.g. \f$\beta_2 \geq 3\f$ and \f$\beta_2 \leq 2\f$)
         * will result in drawRejection() never succeeding, and an error from drawGibbs().
         */
        RestrictionProxy upperBound(size_t k);

        /** Const version of above. */
        const RestrictionProxy upperBound(size_t k) const;

        /** Returns a proxy object that allows assigning coefficient lower bounds.  Element `i` is
         * the lower bound for `beta[i]`, for `i` values between 0 and K-1.
         *
         * For example, to add the restriction \f$\beta_2 \geq 0\f$:
         * 
         *     model.lowerBound(2) = 0.0;
         *
         * Warning: no checking is performed for the viability of specified bounds: specifying
         * conflicting upper and lower bounds (e.g. \f$\beta_2 \geq 3\f$ and \f$\beta_2 \leq 2\f$)
         * will result in drawRejection() never succeeding, and an error from drawGibbs().
         */
        RestrictionProxy lowerBound(size_t k);

        /** Const version of above. */
        const RestrictionProxy lowerBound(size_t k) const;

        /** Returns a proxy object that allows adding both upper and lower bounds for an individual
         * parameter.  Element `i` is the restriction object for `beta[i]`, for `i` values between 0
         * and K-1.
         *
         * To add a restriction, use the <= or >= operator.
         *
         * For example, to add the restrictions \f$0 \leq \beta_2 \leq 5\f$:
         *
         *     model.restrict(2) >= 0;
         *     model.restrict(2) <= 5;
         *
         * The restriction object's <= and >= operators return the object itself, so you can chain
         * restrictions together, such as the following (equivalent to the above):
         *
         *     model.restrict(2) >= 0 <= 5;
         */
        RestrictionIneqProxy restrict(size_t k);

        /** Const version of above. */
        const RestrictionIneqProxy restrict(size_t k) const;

        /** Adds a restriction of the form \f$R\beta \leq r\f$, where \f$R\f$ is a 1-by-`K()` row
         * vector selecting \f$\beta\f$ elements.  For example, to add a \f$\beta_2 \in [-1, 3.5]\f$
         * restriction, any of the following two approaches can be used:
         *
         *     model.upperBound(2) = 3.5;
         *     model.lowerBound(2) = -1;
         *
         *     model.restrict(2) >= -1 <= 3.5;
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

        /** The different draw modes supported by the `draw_mode` parameter and the `draw(mode)`
         * method.
         */
        enum class DrawMode { Auto, Gibbs, Rejection };

        /// The default draw mode used when calling draw() without an argument.
        DrawMode draw_mode = DrawMode::Auto;

        /** The last actual draw mode used, either DrawMode::Gibbs or DrawMode::Rejection.
         * DrawMode::Auto is only set before any draw has been performed.  If DrawMode::Auto was
         * used for the last draw, this value can be used to detect which type of draw actually took
         * place.
         */
        DrawMode last_draw_mode = DrawMode::Auto;

        /** Performs a draw by calling `draw(draw_mode)`.
         *
         * \sa draw(mode)
         * \sa draw_mode
         */
        const Eigen::VectorXd& draw() override;

        /** Perfoms a draw using the requested draw mode, which must be one of:
         *
         *     - DrawMode::Gibbs - for gibbs sampling.
         *     - DrawMode::Rejection - for rejection sampling.
         *     - DrawMode::Auto - try rejection sampling first, switch to Gibbs sampling if required
         *
         * When using DrawMode::Auto, rejection sampling is tried first; if the rejection sampling
         * acceptance rate falls below `draw_auto_min_success_rate` (0.2 by default), rejection
         * sampling is aborted and the current and all subsequent calls to draw() will use Gibbs
         * sampling.  (Note that switching also requires a minimum of `draw_rejection_max_discards`
         * (default: 100) draw attempts before switching).
         */
        const Eigen::VectorXd& draw(DrawMode mode);

        /** Draws a set of parameter values, incorporating the restrictions using Gibbs sampling
         * rather than rejection sampling (as done by drawRejection()).  Returns a const reference
         * to the `k+1` length vector where elements `0` through `k-1` are the \f$\beta\f$ values
         * and element `k` is the \f$h^{-1} = sigma^2\f$ draw.
         *
         * The Gibbs sampler draw uses truncated normal drawing as described in Rodriguez-Yam,
         * Davis, and Scharf (2004), "Efficient Gibbs Sampling of Truncated Multivariate Normal with
         * Application to Constrained Linear Regression," with the additional steps described below.
         *
         * The behaviour of drawGibbs() is affected by two parameters:
         * - `draw_gibbs_burnin` controls the number of burn-in draws to perform for the first
         *   drawGibbs() call
         * - `draw_gibbs_thinning` controls how many (internal) draws to skip between subsequent
         *   drawGibbs() draws.
         *
         * \throws draw_failure if called on a model with no restrictions
         *
         * \par Implementation details: initial value
         *
         * Rodriguez-Yam et al.\ do not discuss the initial draw, but this seems a large omission
         * from the paper: unlike Geweke's earlier work (which is their main benchmark for
         * comparison), their algorithm can fail for points outside the initial truncated region.
         * (Geweke's algorithm, in contrast, does not, though it may wander around outside the
         * truncated region for many initial draws, and may thus require a larger burn-in period).
         * Specifically, in iteration \f$t\f$ for parameter \f$j\f$ in the Rodriguez-Yam et al.\
         * algorithm, one must identify the lower and upper-bounds for parameter \f$j\f$ as required
         * by the model restrictions conditional on parameters \f$1, \ldots, j-1\f$ from the current
         * iteration and \f$j+1, \ldots, k\f$ from the previous iteration.  If there are linear
         * restrictions involving multiple parameters that are violated by the previous draw (or
         * starting point), it can easily be the case that the upper bound on parameter \f$j\f$ is
         * below the lower bound.
         *
         * For example, consider the restrictions \f$\beta_1 \geq 0\f$ and \f$\beta_1 \leq
         * \beta_2\f$, with starting location \f$\beta = [-0.5, -1]^\top\f$.  From the first
         * restriction, we get a lower limit for \f$\beta_1^{(1)}\f$ of 0; applying the second
         * restriction \f$\beta_1^{(1)} \leq \beta_2^{(0)} = -1\f$, in other words, an upper limit
         * of -1.  The algorithm of Rodriguez-Yam et al.\ then requires a draw of
         * \f$\beta_1^{(1)}\f$ truncated to satisfy both \f$\beta_1 \leq -1\f$ and \f$\beta_1 \geq
         * 0\f$, which is clearly impossible.
         *
         * To solve the above problem, we thus need to be certain that the constraints are
         * satisfied.  It doesn't particularly matter how they are satisfied: because of the burn-in
         * period of the sampler, the effects of a initial value far off in the tail of the
         * distribution aren't particularly relevant.  Thus, to choose an initial point (if one
         * hasn't been explicitly given by a call to gibbsInitialize()) we simply take a draw from
         * the unrestricted distribution, then pass it to gibbsInitialize() to manipulate that point
         * into the restriction-satisfying parameter space.  If that fails (specifically, after
         * \f$10k\f$ steps it still hasn't found a suitable point that doesn't violate
         * restrictions), another point is drawn and passed to gibbsInitialize().  This is repeated
         * for a total of 10 attempts, after which we give up with an error (as this most likely
         * indicates that the model has impossible-to-satisfy restrictions).
         *
         * Once an acceptable initial position is found, the algoritm followed is Rodriguez-Yam et
         * al. with one addition described below, subject to the `draw_gibbs_burnin` and
         * `draw_gibbs_thinning` parameters.
         *
         * \par Implementation details: drawing h values
         *
         * Rodriguez-Yam et al.\ is primarily concerned with sampling from a truncated normal
         * distribution, and doesn't satisfactorily describe drawing suitable \f$\sigma^2\f$ values.
         * The approach of Geweke (1991), which addresses drawing from a truncated multivariate
         * \f$\f$-distribution, is essentially to use rejection sampling to take draws from the
         * unrestricted Gamma distribution, accepting the first draw that, had we used it instead of
         * the \f$\sigma^2\f$ from the previous iteration, would not have violated the constraints.
         * Though that approach could be used here (though not exactly as Geweke states it, because
         * we don't require the \f$R\f$ restriction matrix to be invertible, unlike Geweke), it
         * seems preferable to use a deterministic procedure to obtain the next \f$\sigma^2\f$ draw.
         * The procedure is as follows.
         *
         * To incorporate the notation of Rodriguez-Yam et al., let \f$A = L^{-1}\f$, where
         * \f$LL^\top = V\f$, where \f$\overline{s}^2 V\f$ is the (non-truncated) posterior
         * covariance matrix.  For draw \f$t\f$, the algorithm has just drawn a set of values from
         * \f$Z \sim N_S\left(\alpha, \left(\sigma^2\right)^{(t-1)}\right)\f$, where \f$\alpha
         * \equiv A\overline\beta\f$ and \f$S\f$ is the truncation region defined by the model's
         * linear restrictions.  What we need is to find a range for
         * \f$\left(\sigma^2\right)^{(t)}\f$ such that, had we taken exactly the same \f$Z\f$ draw
         * but scaled it to have variance \f$\left(\sigma^2\right)^{(t)}I\f$ instead of
         * \f$\left(\sigma^2\right)^{(t-1)}I\f$, we would not have violated any constraints.
         *
         * To do this, start with \f$l = 0, u=\infty\f$ (i.e. no restrictions), then consider
         * the linear restrictions one at a time.  For restriction \f$i\f$, we have: \f$R_i A^{-1}
         * \left(\alpha + \sigma Y\right) \leq r_i\f$, where \f$Y \equiv \frac{Z - \alpha}{\sigma}\f$, so that
         * the value inside the brackets is just \f$Z\f$, rewritten as the mean plus a
         * scaled draw from a multivariate \f$N(0, I_k)\f$.  Rewriting the constraint, we get:
         * \f$\sigma R_i A^{-1} Y \leq  r_i - R_i A^{-1} \alpha\f$, which immediate leads to either
         * \f$\sigma \leq \frac{r_i - R_i A^{-1} \alpha}{R_i A^{-1} Y}\f$ (if \f$R_i A^{-1} Y\f$ is
         * positive) or \f$\sigma \geq \frac{r_i - R_i A^{-1} \alpha}{R_i A^{-1} Y}\f$ (if \f$R_i
         * A^{-1} Y\f$ is negative).  If this bound is more binding that the current \f$[l, u]\f$
         * bound, update \f$l\f$ or \f$u\f$ appropriately.
         *
         * Finally, convert the \f$[l, u]\f$ bounds on the \f$\sigma\f$ into bounds on \f$h\f$ of
         * \f$[u^{-2}, l^{-2}]\f$, then convert these into the cdf values of the gamma distribution
         * for \f$h\f$, denoting these bounds \f$\alpha\f$ and \f$\omega\f$.  Now we draw
         * \f$\upsilon \sim U[\alpha, \omega]\f$, then pass this through the inverse cdf of the gamma
         * distribution.  This yields \f$h\f$ of the desired distribution, from which
         * \f$\left(\sigma^2\right)^{(t)} = h^{-1}\f$ is our desired \f$\sigma^2\f$ draw.
         */
        virtual const Eigen::VectorXd& drawGibbs();

        /** Initialize the Gibbs sampler with the given initial values of beta, adjusting it, if
         * necessary, to satisfy the model constraints.  If drawGibbs() has previously been called,
         * calling this method will cause the next drawGibbs() call to perform a burn-in.  If the
         * initial value satisfies all the restrictions, it is used as is.  Otherwise, the starting
         * point is adjusted with respect to the violated constraints repeatedly until an admissable
         * point is found, or the maximum number of adjustments is reached.  You can access the
         * final, possibly adjusted value by calling gibbsLast().
         *
         * The specific adjustment is as follows:
         * - Let \f$v \equiv R\beta - r\f$ be the vector of violation amounts, and \f$v_+ \equiv
         *   \sum_{i=1}{l}I(v_i > 0)v_i\f$, i.e.\ the sum of the strictly positive values of
         *   \f$v\f$.
         * - Select a random row (with equal probability of each) from the set of rows of \f$v\f$
         *   that are strictly positive (i.e.  from the violated constraints).
         * - Adjust the current position by moving directly towards the nearest point on the
         *   constraint, and move a distance of 1.5 times the distance to that nearest point.  In
         *   other words, we "overshoot" the constraint boundary by 50%.
         * - Repeat until either all the constraints are satisfied.  If the constraints still
         *   haven't been satisfied after 100 iterations, abort with an exception.
         *
         * Notes:
         * - The reason for the overshoot is to attempt to avoid a condition where two constraints
         *   are at an acute angle to each other; depending on the initial violation (in particular,
         *   in any initial position where jumping to the boundary of either constraint will violate
         *   the other constraint), we could end up "stairstepping" towards the constraint
         *   intersection point, but require (mathematically) an infinite number of steps to reach
         *   it; with numerical imprecision we might or might not reach it at all.  The overshoot is
         *   designed to help alleviate this.
         * - This algorithm is invariant to scaling in the restriction specification (that is,
         *   restriction specifications of \f$\beta_1 + \beta_2 \leq 2\f$ and \f$5\beta_1 + 5\beta_2
         *   \leq 10\f$ are treated identically).
         * - This algorithm specifically does not select violated rows based on the size of the
         *   violation or the distance from current position to the violated boundary, because there
         *   is no particular reason to think that cross-row restrictions are comparable.  For
         *   example, suppose a model has restriction \f$\beta_1 \geq 1\f$ with data \f$X_1\f$ and
         *   another model uses the rescaled data \f$10X_1\f$ with restriction \f$\gamma_1 \geq
         *   10\f$.  The size of the violation from the same (relative) draw, and the Euclidean
         *   distance from the inadmissable point to the restriction boundary, will have been scaled
         *   by 10; thus depending on either in the selection of constraints to attempt to satisfy
         *   would result in model reparameterizations fundamentally changing the algorithm, which
         *   seems an undesirable characteristic.
         * - One potential solution to the above would be to normalize the parameters.  This,
         *   however, would require that the number of parameters involved in restrictions equals
         *   the number of restrictions (so that the restricted values can be normalized by the same
         *   matrix that normalizes the restricted \f$\beta\f$'s), but that sort of limitation isn't
         *   otherwise required by the Gibbs algorithm used in this class.
         *
         * \param initial the initial value of \f$\beta\f$.  May also include an extra
         * \f$\sigma^2\f$ value, which will be ignored.
         * \param max_tries the maximum tries to massage the initial point into unrestricted
         * parameter space.  Defaults to 100.
         * \throws constraint_failure if the constraint cannot be satisfied after `max_tries`
         * adjustments.
         * \throws logic_error if called with a vector with size != K
         */
        void gibbsInitialize(const Eigen::Ref<const Eigen::VectorXd> &initial, unsigned long max_tries = 100);

        /** Returns the beta values associated with the current Gibbs sampler value.  Note that this
         * is not necessarily a draw: in particular, the current value might be an adjusted initial
         * value.  For anything other than obtaining the initial value, calling lastDraw() is a
         * better choice (particularly because it doesn't require any calculations or vector copy,
         * unlike this method).
         *
         * If there is no last value at all, this returns an empty vector.
         */
        Eigen::VectorXd gibbsLast();

        /** Exception class used to indicate that one or more constraints couldn't be satisfied.
         */
        class constraint_failure : public draw_failure {
            public:
                /** Constructor.
                 * \param what the exception message.
                 */
                constraint_failure(const std::string &what) : draw_failure(what) {}
        };

        /** Returns a draw from a truncated univariate distribution given the truncation points, the
         * cdf and inverse cdf functions, cdf complement and its inverse cdf complement, and
         * (optionally) the median.  The complements are used to ensure better numerical precision
         * for regions of the distribution with cdf values above 0.5.
         *
         * Note: it is recommended to check that the truncation bounds of the distribution are
         * actually within the support of the distribution; if neither is, it is more efficient and
         * more accurate to sample from the non-truncated distribution instead.  It may also be more
         * efficient (probabilistically) to use rejection sampling when the truncation bounds are
         * reasonably far out in opposite tails of the distribution.
         *
         * \param min the truncated region lower bound
         * \param max the truncated region upper bound
         * \param cdf returns the cdf value for a distribution value.
         * \param quantile returns the distribution value for a cdf value (i.e. the inverse cdf)
         * \param cdf_complement returns 1 minus the cdf, but typically with better precision.
         * \param quantile_complement returns the quantile from a cdf complement value.
         * \param median the median of the distribution, above which cdf complements will be used
         * instead of cdf values, for increased precision.  If omitted, quantile will be used first,
         * and, if the result is greater than 0.5, quantile_complement called instead.  (This
         * double-lookup is done at most once, either for min or max but not both).  If the median
         * can be calculated without an inverse cdf lookup, passing it is recommended (and strongly
         * recommended when inverse cdf lookups are expensive) as it often avoids an unnecessary
         * inverse cdf (or inverse complement cdf) calculation whenever the median is not in the
         * truncated region.
         */
        static double truncDist(
                double min,
                double max,
                const std::function<double(double)> &cdf,
                const std::function<double(double)> &cdf_complement,
                const std::function<double(double)> &quantile,
                const std::function<double(double)> &quantile_complement,
                double median = std::numeric_limits<double>::signaling_NaN());

        /** Draws a truncated normal with mean parameter `mean`, standard deviation `sd`, truncated
         * to the range `[min,max]`.  This is done by drawing from a
         * \f$U[\Phi_{mean,sd}(0),\Phi_{mean,sd}(1)]\f$ where \f$Phi_{mean,sd}\f$ is the cdf of a
         * normal with mean $mean$ and standard deviation $sd$.  The resulting uniform value is then
         * passed through the normal cdf inverse (using boost's erfc_inv()) to get the truncated
         * draw.
         *
         * If min is negative infinity and max is positive infinity, the truncation is bypassed
         * entirely and simply `mean + sd * eris::Random::rstdnorm()` is returned.
         *
         * \throws constraint_failure if `min >= max`
         */
        static double truncNorm(double mean, double sd, double min, double max);

        /** Draws a truncated gamma with shape parameter `shape` and scale parameter `scale`,
         * truncated to the range `[min,max]`.
         *
         * This is done just like truncNorm, but using the cdf and inverse cdf for a gamma
         * distribution instead of normal distribution.
         *
         * If `min` is 0 (or negative) and `max` is infinity, this simply returns a non-truncated
         * gamma distribution draw, bypassing the above cdf/inverse cdf calculations.
         *
         * \throws constraint_failure if `min >= max` or `max <= 0`
         */
        static double truncGamma(double shape, double scale, double min, double max);

        /** This method attempts to draw the `K+1` \f$(\beta,h^{-1})\f$ values as if no restrictions
         * are in place (by calling Linear::draw()), then repeats the draw as long as the drawn
         * values violate one or more of the models restrictions (or `draw_discards_max`
         * unsuccessful draws occur).
         *
         * This is more efficient than drawGibbs() when the likelihood of restricted values being
         * drawn is very small; it is usually recommended to call `draw()` instead, which tries this
         * and switches to Gibbs if the acceptance rate is low.
         *
         * After calling, whether returning or throwing a draw_failure exception,
         * `draw_rejection_discards`, `draw_rejection_discards_cumulative`, and
         * `draw_rejection_success_cumulative` will have been updated to reflect the draw
         * statistics.
         *
         * \param max_discards if specified and positive, the maximum number of failed draws before
         * aborting (by throwing an exception).  The default value, -1, uses
         * `draw_rejection_max_discards`.
         *
         * \throws Linear::draw_failure exception if, due to the imposed restrictions, at least
         * `draw_discards_max` sequential draws are discarded without finding a single admissible
         * draw.
         */
        virtual const Eigen::VectorXd& drawRejection(long max_discards = -1);

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

        int draw_rejection_discards_last = 0, ///< Tracks the number of inadmissable draws by the most recent call to drawRejection()
            draw_rejection_success = 0, ///< The cumulative number of successful rejection draws
            draw_rejection_discards = 0, ///< The cumulative number of inadmissable rejection draws
            draw_rejection_max_discards = 100, ///< The maximum number of inadmissable draws for a single rejection draw before aborting
            draw_gibbs_burnin = 100, ///< The number of burn-in draws for the first Gibbs sampler draw
            draw_gibbs_thinning = 3; ///< drawGibbs() uses every `draw_gibbs_thinning`th sample (1 = use every draw)
        double draw_auto_min_success_rate = 0.2; ///< The minimum draw success rate below which we switch to Gibbs sampling

        /// Accesses the restriction coefficient selection matrix (the \f$R\f$ in \f$R\beta <= r\f$).
        Eigen::Ref<const Eigen::MatrixXd> R() const;

        /// Accesses the restriction value vector (the \f$r\f$ in \f$R\beta <= r\f$).
        Eigen::Ref<const Eigen::VectorXd> r() const;

        /** Overloaded to append the restrictions after the regular Linear details.
         */
        virtual operator std::string() const override;

        /** The display name of the model to use when printing it.  Defaults to "Linear" but
         * subclasses should override.
         */
        virtual std::string display_name() const;

    protected:
        /// Creates a LinearRestricted from a Linear rvalue
        LinearRestricted(Linear &&move) : Linear(std::move(move)) {}

        /** Resets any drawn values and draw-related variables.  Called automatically when adding a
         * restriction, in addition to the `Linear` cases (updating and weakening).
         */
        virtual void reset() override;

        /** Overridden to also reset mean_beta_draws_ to 0 (to force predict() value redrawing). */
        virtual void discardForce(unsigned int burn) override;

        /** Proxy object that converts assignments into restriction rows on the associated
         * LinearRestricted model.
         *
         * \sa upperBound()
         * \sa lowerBound()
         */
        class RestrictionProxy final {
            public:
                /** Adds a restriction on the referenced parameter.
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
                friend class LinearRestricted;
                LinearRestricted &lr_;
                const size_t k_;
                const bool upper_;
        };

        /** Proxy object that converts <= and >= inequality operators into new restriction rows on
         * the associated LinearRestricted model.
         *
         * \sa restrict()
         */
        class RestrictionIneqProxy final {
            public:
                /** Adds an upper-bound restriction on the referenced parameter.
                 */
                RestrictionIneqProxy& operator<=(double r);
                /** Adds a lower-bound restriction on the referenced parameter.
                 */
                RestrictionIneqProxy& operator>=(double r);
                /** Returns true if there is an upper-bound restriction on the associated parameter.
                 * This method works whether the restriction was added via a proxy object or as a
                 * single-parameter restriction given to addRestriction() (or variants).
                 *
                 * Specifically this looks for any restriction with a positive selector for the
                 * requested coefficient and 0 for all other coefficient selectors.
                 */
                bool hasUpperBound() const;
                /** Returns the most-binding upper-bound restriction on the associated parameter.
                 * If the parameter has no upper-bound restrictions at all, returns NaN.
                 */
                double upperBound() const;
                /** Returns true if there is an lower-bound restriction on the associated parameter.
                 * This method works whether the restriction was added via a proxy object or as a
                 * single-parameter restriction given to addRestriction() (or variants).
                 *
                 * Specifically this looks for any restriction with a negative selector for the
                 * requested coefficient and 0 for all other coefficient selectors.
                 */
                bool hasLowerBound() const;
                /** Returns the most-binding lower-bound restriction on the associated parameter.
                 * If the parameter has no lower-bound restrictions at all, returns NaN.
                 */
                double lowerBound() const;
            private:
                /** Creates a new RestrictionProxy object that, when assigned to, adds an upper or
                 * lower-bound restriction on beta[`k`].
                 *
                 * \param the LinearRestricted object where the restriction will be added;
                 * \param k the parameter index
                 */
                RestrictionIneqProxy(LinearRestricted &lr, size_t k);
                friend class LinearRestricted;
                LinearRestricted &lr_;
                const size_t k_;
        };

        /** Returns true if there is a single-parameter, upper- or lower-bound restriction on
         * `beta[k]`.  Note that this ignores any multi-parameter restriction involving `beta[k]`.
         *
         * \sa RestrictionProxy::restricted()
         */
        bool hasRestriction(size_t k, bool upper) const;

        /** Returns the most-binding single-parameter, upper- or lower-bound restriction on
         * `beta[k]`, or NaN if `beta[k]` has no such restriction.  Note that this ignores any
         * multi-parameter restriction involving `beta[k]`.
         *
         * \sa RestrictionProxy::operator double()
         */
        double getRestriction(size_t k, bool upper) const;

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

        /** The cache of drawn beta vectors used for prediction.  Must not be used if
         * mean_beta_draws_ is 0.
         */
        Eigen::VectorXd mean_beta_;

        /// The number of beta draws used to calculate mean_beta_
        long mean_beta_draws_ = 0;

    private:
        // Values used for Gibbs sampling.  These aren't set until first needed.
        std::shared_ptr<Eigen::MatrixXd> gibbs_D_; // D = R A^{-1}
        std::shared_ptr<Eigen::VectorXd> gibbs_alpha_, gibbs_last_; // alpha = A \mu; last = Z (not X!)
        long gibbs_draws_ = 0;
};

}}
