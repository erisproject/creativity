#pragma once
#include <cstdint>

namespace creativity {

/** Simulation parameters that are used to configure the simulation when calling setup(). */
struct CreativitySettings {
    /// The number of readers in the simulation
    uint32_t readers = 100;

    /// The number of dimensions in the simulation
    uint32_t dimensions = 2;

    /** The density of the simulation in average number of readers per \f$unit^{dimensions}\f$.
     * This, combined with `readers` and `dimensions` implicitly defines the boundaries where the
     * simulation space wraps to the opposite boundary.
     *
     * For example, for `readers = 100, dimensions = 2, density = 1`, the boundaries will be at ±5
     * in each dimensions.  `readers = 100, dimensions = 3, density = 2` would result in boundaries
     * at ±1.842 in each of the three dimensions.
     *
     * This field is used instead of boundary if `use_density` is true (the default).
     */
    double density = 1.0;

    /** The boundary in each dimension to use.  This is used instead of `density` if `use_density`
     * is false; if it is true, this is updated to the boundary implied by the `density` setting
     * when the simulation begins.
     */
    double boundary = 5.0;

    /** When setting up the simulation, if this is true (the default), use `density` to calculate
     * the boundary (and then replace `boundary` with the calculated value).  If false, use
     * `boundary` directly, and replace `density` with the appropriately calculated value.
     *
     * Once a simulation is configured (or loaded from a file) both density and boundary will agree,
     * and so this parameter is irrelevant.
     */
    bool use_density = true;

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
     *
     * Like `cost_fixed` and `cost_unit` this only specifies the default reader income: reader
     * income can be updated on a individual reader level.
     */
    double income = 1000.0;

    /** The period in which the sharing network is introduced.  If 0, sharing is never invented.
     * (Since t=0 is the initial setup, without any actions, set to 1 to having sharing available
     * immediately).
     */
    uint64_t sharing_begins = 100;

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
    double sharing_link_proportion = 0.1;

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
        /** Initial book prices are distributed Uniform[`p_min`, `p_max`] */
        double p_min = 2.0;
        /** Initial book prices are distributed Uniform[`p_min`, `p_max`] */
        double p_max = 5.0;
        /** The probability of a book being kept on the market for another period. */
        double prob_keep = 0.5;
        /** The relative price-above-marginal-cost of a book being kept on the market for another
         * period.  For example, if this is 0.75, the current price is 3, and `parameters.cost_unit`
         * is 1, then the new price will be 1 + 0.75*(3-1) = 2.5.
         */
        double keep_price = 0.5;
    } initial;
};

double boundaryFromDensity(uint32_t readers, uint32_t dimensions, double density);
double densityFromBoundary(uint32_t readers, uint32_t dimensions, double boundary);

}
