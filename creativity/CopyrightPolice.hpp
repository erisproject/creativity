#pragma once
#include "creativity/Creativity.hpp"
#include <eris/Agent.hpp>
#include <eris/Optimize.hpp>
#include <boost/math/distributions/normal.hpp>

namespace creativity {

class Creativity;

/** This class represents a copyright policing agent that (probabilistically) catches readers who
 * obtain books via piracy, fining them and redistributing fines to infringed-upon authors.
 */
class CopyrightPolice : public eris::Agent, public virtual eris::interopt::Apply, public virtual eris::intraopt::Finish {
public:
    CopyrightPolice() = delete; ///< Not default constructible
    /** Constructor for a new CopyrightPolice agent.  Takes a reference to the creativity object
     * and the lump sum tax amount.
     *
     * \throws std::domain_error if catch tax is negative, or if sigma is negative for the level
     * of the catch tax.
     */
    explicit CopyrightPolice(const Creativity &creativity);

    /// When the period advances, we take the lump sum tax from all agents.
    void interApply() override;

    /// Override priority to run after the Reader's interApply has deposited income.
    double interApplyPriority() const override { return 1.0; }

    /** When the period finishes, we probabilistically detect piracy; readers found pirating and
     * fined, and the fines redistributed to the authors of the pirated works.
     */
    void intraFinish() override;

    /// Returns the lump sum per-reader tax collected each period
    double tax() const { return creativity_.parameters.policy_catch_tax; }

    /** Returns the \f$\mu\f$ parameter; the probability of being caught is the CDF of a
     * \f$\mathcal{N}(\mu, \sigma^2)\f$ distribution.
     */
    double mu() const { return normal_.mean(); }

    /** Returns the \f$\sigma\f$ parameter; the probability of being caught is the CDF of a
     * \f$\mathcal{N}(\mu, \sigma^2)\f$ distribution.
     */
    double sigma() const { return normal_.standard_deviation(); }

    /** Returns the probability of being caught for a given number of pirated books.  This is
     * simply the CDF of a normal distribution with parameters as given by mu() and sigma().
     */
    double prob(unsigned pirated) const { return cdf(normal_, pirated); }

    /** Returns the cost of being accused of piracy (whether guilty or not).  This is simply a
     * shortcut for `creativity_.parameters.policy_catch_cost`.
     *
     * Unlike fine, this isn't a paid amount, but rather represents the opportunity cost to the user
     * of being caught (as such it isn't really a CopyrightPolice property, but is provided here
     * nonetheless for convenience).
     */
    double cost() const { return creativity_.parameters.policy_catch_cost; }

    /** Returns the fine for being caught pirating `pirated` books this period.  Will always return
     * 0 for 0 pirated books, otherwise evaluates the `creativity_.parameters.policy_catch_fine`
     * polynomial to determine the fine.  Returns 0 if the polynomial evaluates to a negative
     * number.
     */
    double fine(unsigned pirated);

protected:
    /// The Creativity object that owns the simulation this reader belongs to
    const Creativity &creativity_;

    /// The lump sum tax amount
    double tax_;

    /// The normal distribution used to generate probabilities
    boost::math::normal normal_;
};

}
