#pragma once
#include <eris/types.hpp>
#include <cstdint>
#include <array>

namespace creativity {

/** Simulation parameters that are used to configure the simulation when calling setup(). */
struct CreativitySettings {
    /// The number of readers in the simulation
    uint32_t readers = 100;

    /// The number of dimensions in the simulation
    uint32_t dimensions = 2;

    /** The boundary in each dimension to use for simulation members.  Must be a positive value.
     * The the exact boundary is calculated from the simulation density at the initialization of the
     * simulation as the density could theoretically change during the simulation, but the boundary
     * defines the world and so is fixed forever.
     */
    double boundary = 5.0;

    /** An authored book will be located at a distance of \f$mx\f$ from the author's location at the
     * time of writing, where this value is \f$m\f$ and \f$x \sim \chi^2_1\f$.  This distance thus
     * has the given mean and a standard deviation of \f$\sqrt{2}m\f$.
     *
     * The final location is drawn from a uniform distribution over all points with the given
     * distance from the author.
     *
     * If this value is 0, the book is located exactly at the author's position at the time of
     * writing and no draw is performed.
     */
    double book_distance_mean = 0.0;

    /** The standard deviation of a book quality draw.  When a reader obtains a book, his subjective
     * quality is drawn from \f$N(Q, s)\f$, where \f$Q\f$ is the book's base quality as decided by
     * the author and \f$s\f$ is this setting.
     *
     * As long as this setting is positive, this means book quality is subjective; if 0, all readers
     * perceive the book as having the same quality.
     */
    double book_quality_sd = 1.0;

    /** Between periods, readers take a random step of length \f$mx\f$, where this value is \f$m\f$
     * and \f$x \sim \chi^2_1\f$.  This value has mean \f$m\f$ and standard deviation \f$\sqrt{2}
     * m\f$.  If 0, readers remain at their initial positions forever.
     *
     * The final location is selected from a uniform distribution over all points with the given
     * distance from the reader's initial position.
     *
     * The default is 0.5.
     */
    double reader_step_mean = 0.5;

    /** The "shape" parameter \f$\beta\f$ of the reader effort.
     *
     * \sa Reader.creation_shape
     */
    double reader_creation_shape = 0;

    /** Reader.creation_scale values are drawn from \f$U[a,b]\f$, where this value is the lower
     * bound.  Defaults to 0.
     *
     * \sa Reader.creation_scale
     */
    double reader_creation_scale_min = 0.0;

    /** Reader.creation_scale values are drawn from \f$U[a,b]\f$, where this value is the size of
     * the range `b-a`.  Defaults to 10 (so the values are drawn from U[0,10]).
     */
    double reader_creation_scale_range = 10.0;

    /** The average length of time (in simulation periods) it takes to create a book.  If 0, books
     * are created instantly; otherwise books take this value -1/+0/+1 with equal probability of
     * each (so that this is still the mean creation time).  (The randomness helps avoid strongly
     * cyclical authorship behaviour)
     */
    uint32_t creation_time = 3;

    /** The fixed cost component of initial creation.  This plus the expended effort level is
     * incurred to create.
     */
    double creation_fixed = 100.0;

    /** The fixed cost of keeping a book on the market.
     */
    double cost_market = 10.0;

    /** The unit cost of selling a copy of a book (note that only copies of books still on the
     * market can be sold), which is also the cost of obtaining a pirated copy (but is then incurred
     * by the recipient, rather than the author).
     *
     * Note that, if changing this, it is advisable to also consider changing
     * `parameters.initial.p_min`.
     */
    double cost_unit = 1.0;

    /** The cost of receiving a pirated copy of a book, incurred by the recipient.
     */
    double cost_piracy = 1.0;

    /** The per-period external income readers receive.  Effort spent creating a book in a period is
     * subtracted from this amount.
     */
    double income = 1000.0;

    /** The period in which piracy becomes available.  If 0, piracy is never invented.  (Since t=0
     * is the initial setup (without any actions), set to 1 to having sharing become available
     * immediately).
     */
    eris::eris_time_t piracy_begins = 101;

    /** Bitfield specifying the enabled policy response(s) used to address piracy.  Currently
     * supported bit values are:
     * 
     * - 0x1: public tax for public provisioning and sharing; when the policy response begins, a
     *   PublicTracker agent is created that distributes copies of books to anyone at marginal cost.
     *   Authors are compensated in proportion to the number of downloads of their works from a
     *   per-read lump sum tax.
     *
     * - 0x2: public tax for catching people downloading illegal copies.  A lump sum tax is
     *   collected and goes towards the policing effort, with a larger tax resulting in a higher
     *   probability of detection.
     * 
     * Other values are reserved for future use.
     */
    uint32_t policy = 0;

    /** The period in which the policy response to piracy begins.  Which response is used depends on
     * the `policy` parameter.
     */
    eris::eris_time_t policy_begins = 201;

    /** This is the lump-sum tax extracted from each user at the beginning of the period to fund the
     * public sharing mechanism of author compensation.  Has no effect when public sharing is not
     * enabled.
     */
    double policy_public_sharing_tax = 10.0;

    /** This is the lump-sum tax extracted from each user to fund the policing policy to catch and
     * fine readers using piracy.  Has no effect when the catching policy is not enabled.
     */
    double policy_catch_tax = 10.0;

    /** The utility cost of being accused of piracy.  This is incurred by everyone who is accused of
     * piracy, regardless of whether or not they are guilty (since the probability of being caught
     * even with no piracy is > 0).  If guilty, the reader pays an additional fine, determined by
     * `policy_catch_fine`.  This component of the cost is simply a cost, it is not a fine being
     * transferred, unlike `policy_catch_fine`.
     */
    double policy_catch_cost = 10.0;

    /** A person who gets caught and is guilty of piracy pays a fine of \f$c_0 + c_1 #P + c_2
     * #P^2\f$, where \f$#P\f$ is the number of pirated books the person obtained in the caught
     * period.  Thus a caught individual who obtained 3 pirated copies pays \f$c_0 + 3 c_1 + 9
     * c_2\f$.
     *
     * The fine amount is removed from the reader's utility (as a proxy for the reader's wealth,
     * which is outside this model), and is redistributed to the authors of the pirated works.
     * Specifically, the fine is divided into \f$#P\f$ equal shares and each share is transferred to
     * the author of each pirated work.  (It is possible for a single author to receive multiple
     * shares, if multiple copies of his works were pirated by the fined reader).
     *
     * The default is a fine of 50 for the first, +50 for the second, +75 for the third, etc. (so a
     * total fine of 50, 100, 175, 275, ...).
     */
    std::array<double, 3> policy_catch_fine{{25., 12.5, 12.5}};

    /** For the "catch pirates" policy, the probability of being accused of piracy for obtaining
     * \f$x\f$ copies over the network equals the CDF of a \f$\mathcal{N}(\mu, \sigma^2)\f$
     * distribution evaluated at \f$x\f$; that is, \f$\Phi\(\frac{x - \mu}{\sigma}\)\f$.
     *
     * \f$\mu\f$ and \f$\sigma\f$ themselves are quadratic functions of the per-user piracy
     * detection expenditure, i.e., functions of `policy_response_tax`.
     *
     * This array contains the constant, linear, and quadratic coefficients that define \f$\mu = c_0
     * + c_1 \tau + c_2 \tau^2\f$, where \f$\tau\f$ is the `policy_response_tax` value.  Note that
     * the probability of detection equals 50% at \f$\mu\f$, 5% at \f$\mu - 1.96 \sigma\f$, 95%
     * at \f$\mu + 1.96\f$.  The default is \f$\mu = 10 - 0.5 \tau\f$.
     *
     * Note that the probability of being caught is never 0, even for an individual who does not
     * pirate; such false accusations can still incur a cost, or can be ignored, depending on the
     * `policy_catch_cost` variable.
     *
     * The default for this an mu are selected with a range of tau of (0, 100] in mind: close to 0,
     * the probability of catching someone who doesn't pirate is very close to 0, the 50% detection
     * is at 10 pirated copies in a period, and the 95% is at 16 (i.e. the probability is the CDF of
     * N(10,3)).  At the upper end, the CDF is N(2, 1), i.e. non-pirates are false accussed just
     * under 5% of the time, people pirating 2 are caught 50% of the time, and 95% is hit after
     * pirating 4 copies.
     */
    std::array<double, 3> policy_catch_mu{{10., -.08, 0.}};

    /** Like `policy_catch_mu`, but for the quadratic coefficients of \f$\sigma = c_0 +
     * c_1 \tau + c_2 \tau^2\f$.  The default is \f$\sigma = 3 - .02 \tau\f$.  Note that this will
     * fail (immediately) if tau is >= 150.
     *
     * \see policy_catch_mu
     */
    std::array<double, 3> policy_catch_sigma{{3., -.02, 0.}};

    /** The number of sharing/friendship links as a proportion of the maxinum number of sharing
     * links possible (which is \f$\frac{R(R-1)}{2}\f$, where \f$R\f$ is the number of readers).
     *
     * Links as assigned randomly between agents when the simulation is initially set up, with an
     * equal probably of each possible link being selected.
     *
     * The default is 10% coverage of maximum potential links (rounded to the nearest integer).  In
     * the default 100-reader simulation, this is 495 links.
     *
     * The value must be in \f$[0, 1]\f$.
     */
    double piracy_link_proportion = 0.1;

    /** The factor by which parameter standard deviations should be multiplied when using it as a
     * prior for a subsequent period, after the initial burn-in period.
     */
    double prior_scale = 1.02;

    /** The factor by which to multiply standard deviations in the first piracy period.  This value
     * overrides `prior_scale` in the `piracy_begins` period.
     */
    double prior_scale_piracy = 2;

    /** The factor by which to multiply standard deviations in the first policy response period.
     * This value overrides `prior_scale` in the `policy_response_begins` period.
     */
    double prior_scale_policy = 2;

    /** The prior weight to use during the burn-in period; typically much larger than prior_scale
     * so that the simulation results of pre-belief initial parameters have significantly less
     * weight in readers' beliefs.
     */
    double prior_scale_burnin = 1.5;

    /** The number of "burn-in" periods, during which priors are discounted at a higher rate. */
    uint32_t burnin_periods = 20;

    /** The number of model draws to use for prediction.  Higher values yield more "accurate"
     * predictions, but lower values may be desirable to introduce more random agent behaviour.
     */
    uint32_t prediction_draws = 100;

    /** The values in this struct define fixed probabilities and distributions of simulation
     * actions.  This is needed because, in the initial simulation periods, readers only have
     * noninformative beliefs which are completely useless for predicting anything in the model.
     *
     * Instead, readers act randomly until useful beliefs are established, according to the values
     * that follow.
     */
    struct {
        /** The probability that a reader writes a book. */
        double prob_write = 0.1;
        /** Initial book creation effort levels are distributed Uniform[`l_min`, `l_min+l_range`].
         * Readers without sufficiently informed quality beliefs use the mean of the resulting
         * quality distribution for initial purchasing decisions. */
        double l_min = 0;
        /** Initial book creation effort levels are distributed Uniform[`l_min`, `l_min+l_range`].
         * Readers without sufficiently informed quality beliefs use the mean of the resulting
         * quality distribution for initial purchasing decisions. */
        double l_range = 100.0;
        /** Initial book prices are distributed `cost_unit +` Uniform[`p_min`, `p_min+p_range`].
         * Note that initial book quality, at default creation shape/scale and l_min/l_range
         * settings, can take a value anywhere from 0 to 46.15, and has a mean of 20.14. */
        double p_min = 4.0;
        /** Initial book prices are distributed `cost_unit +` Uniform[`p_min`, `p_min+p_range`] */
        double p_range = 30.0;
        /** The probability of a book being kept on the market for another period. */
        double prob_keep = 0.5;
        /** The relative price-above-marginal-cost of a book being kept on the market for another
         * period.  For example, if this is 0.75, the current price is 3, and `parameters.cost_unit`
         * is 1, then the new price will be 1 + 0.75*(3-1) = 2.5.
         */
        double keep_price = 0.5;
        /** The minimum number of observations (relative to k) before a belief is used.  Before reaching
         * this threshold, actions are governed by the `initial` values, above.
         *
         * This setting also governs the use of the ProfitStream beliefs of different lengths: authors
         * will not use ProfitStream belief models that do not have at least this number of
         * observations.  (For example, if this setting equals 2 and there are ProfitStream beliefs for
         * 1-, 2-, and 4-period old books, with 15, 7, and 5 observations, only the 1- and 2-period
         * models will be used.)
         *
         * Note that partially informative models (i.e. those with fewer than K linearly-independent
         * observations) are never used, so setting this to a value smaller than 0 has no effect,
         * and often the effective minimum threshold is noticeably higher than 0 (particularly when
         * the model contains dummy or other small integer values).
         */
        int32_t belief_threshold = 5;
    } initial;

};

}

