#pragma once
#include <eris/types.hpp>
#include <cstdint>

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

    /** The standard deviation of a book.  An authored book will be located at a distance drawn from
     * \f$\left|N(0,s)\right|\f$, where \f$s\f$ is this value, in a random direction from the
     * author's location at the time of writing.
     *
     * The direction is drawn from a uniform distribution over the surface of the hypersphere
     * centred on the author with the randomly drawn radius.
     *
     * If this value is 0, the book is located exactly at the author's position at the time of
     * writing.
     */
    double book_distance_sd = 0.5;

    /** The standard deviation of a book quality draw.  When a reader obtains a book, his subjective
     * quality is drawn from \f$N(Q, s)\f$, where \f$Q\f$ is the book's base quality as decided by
     * the author and \f$s\f$ is this setting.
     *
     * As long as this setting is positive, this means book quality is subjective; if 0, all readers
     * perceive the book as having the same quality.
     */
    double book_quality_sd = 1.0;

    /** Between periods, readers take a random step of length \f$N(0, s^2)\f$.  This value is
     * \f$s\f$ (the standard deviation, not the variance \f$s^2\f$).  If 0, readers remain at their
     * initial positions forever.
     *
     * More specifically, a random direction vector of length 1 is multiplied by a draw from this
     * distribution, so negative draws simply reverse the direction of the random vector; hence the
     * actual distance moved is the absolute value of the draw, in a random direction.
     *
     * The default is 0.25.
     */
    double reader_step_sd = 0.25;

    /** The "shape" parameter \f$\beta\f$ of the reader effort.
     *
     * \sa Reader.creation_shape
     */
    double reader_creation_shape = 0;

    /** Reader.creation_scale values are drawn from \f$U[a,b]\f$, where this value is the lower
     * bound.  Defaults to 5.
     *
     * \sa Reader.creation_scale
     */
    double reader_creation_scale_min = 1.0;

    /** Reader.creation_scale values are drawn from \f$U[a,b]\f$, where this value is the upper
     * bound.  Defaults to 15.
     */
    double reader_creation_scale_max = 10.0;

    /** The length of time (in simulation periods) it takes to create a book.  If 0, books are
     * created instantly; if larger, the given number of periods go by before the book is finished.
     */
    uint32_t creation_time = 3;

    /** The fixed cost of keeping a book on the market.
     */
    double cost_fixed = 10.0;

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

    /** The period in which the PublicTracker agent is created.  The public tracker provides
     * marginal cost access to copies of books but taxes all agents a lump sum amount,
     * redistributing the collected tax money to authors proportionally to the number of copies
     * obtained of each author's works.
     */
    eris::eris_time_t public_sharing_begins = 201;

    /** The lump size tax extracted by the public tracker from each agent in each period.
     *
     * Currently this is fixed, but future versions of this code may interpret this value as a
     * starting value and attempt to adjust it to maximize long-term utility.
     */
    double public_sharing_tax = 10.0;

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

    /** The factor by which to multiply standard deviations in the first public sharing period.
     * This value overrides `prior_scale` in the `public_sharing_begins` period.
     */
    double prior_scale_public_sharing = 2;

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
        double prob_write = 0.2;
        /** Initial book qualities are distributed Uniform[`q_min`, `q_max`].  Uninformed readers
         * draw from this same distribution when guessing quality. */
        double q_min = 0;
        /** Initial book qualities are distributed Uniform[`q_min`, `q_max`].  Uninformed readers
         * draw from this same distribution when guessing quality. */
        double q_max = 10.0;
        /** Initial book prices are distributed `cost_unit +` Uniform[`p_min`, `p_max`] */
        double p_min = 0.5;
        /** Initial book prices are distributed `cost_unit +` Uniform[`p_min`, `p_max`] */
        double p_max = 5.0;
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

