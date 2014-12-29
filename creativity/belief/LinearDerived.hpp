#pragma once
#include "creativity/belief/Linear.hpp"

namespace creativity { namespace belief {

/** Wrapper around Linear with common methods needed by the Linear-derived classes.  In particular,
 * this handles subclass versions of update() and weaken() that wrap Linear::update and
 * Linear::weaken but return the derived type instead of a Linear.
 *
 * \param D the derived class
 * \param B the base class that D inherits from, which must be either Linear (the default) or an
 * intermediate base class itself derived from Linear (such as LinearRestricted).
 */
template <class D, class B = Linear, typename = typename  std::enable_if<std::is_base_of<Linear, B>::value>::type>
class LinearDerived : public B {
    public:
        /// Inherit all constructors from the base class
        using B::B;

        /** Uses the current object's as a prior to generate a new object whose parameters are the
         * posteriors of this object after incorporating new data.
         *
         * \param y a vector of new y data
         * \param X a matrix of new X data
         */
        [[gnu::warn_unused_result]]
        D update(const Eigen::Ref<const Eigen::VectorXd> &y, const Eigen::Ref<const Eigen::MatrixXd> &X) const & {
            return newDerived(B::update(y, X));
        }

        /** Updates the current object, using its current values as a prior, incorporating the
         * passed-in new data.  This method is called when update() is invoked on an rvalue
         * reference, avoiding requiring an extra copy constructor invocation when chaining
         * operations such as:
         *
         *     model = model.weaken(scale).update(y, X);
         *
         * \param y a vector of new y data
         * \param X a matrix of new X data
         */
        [[gnu::warn_unused_result]]
        D update(const Eigen::Ref<const Eigen::VectorXd> &y, const Eigen::Ref<const Eigen::MatrixXd> &X) && {
            return newDerived(std::move(*this).B::update(y, X));
        }

        /** Weakens a Derived belief, returning a new, weakened Derived object.
         *
         * \sa Linear::weaken
         */
        [[gnu::warn_unused_result]]
        D weaken(double prior_weight) const & {
            return newDerived(B::weaken(prior_weight));
        }

        /** Weakens a Derived belief, updating the current rvalue-reference then returning it.
         *
         * \sa Linear::weaken
         */
        [[gnu::warn_unused_result]]
        D weaken(double prior_weight) && {
            return newDerived(std::move(*this).B::weaken(prior_weight));
        }

    protected:
        /** Called to construct a new Derived object from a Linear base object rvalue.  There is no
         * default implementation.
         */
        virtual D newDerived(Linear &&new_base) const = 0;

        /** Alias for subclasses to refer to this as Parent instead of repeating the entire class
         * with template arguments again.
         */
        using Parent = LinearDerived<D, B>;

        /** Constructor using Linear r-value */
        LinearDerived(Linear &&move) : B(std::move(move)) {}

        /** Default constructor */
        LinearDerived() = default;
};

}}
