#pragma once
#include <list>
#include <memory>
#include <functional>
#include <utility>
#include <ostream>
#include "creativity/data/Variable.hpp"

namespace creativity { namespace data {

/** Class to store a equation.  The class is constructed with the dependent variable, then has an
 * equation added to it using the % operator.  For example:
 *
 *     OLS model(depvar);
 *     model % 1 + var1 + var2;
 *     model % var3;
 *
 * NB: The left-shift operator (<<) would make more logical sense than %, but has the unfortunately
 * pitfall of being of lower precedence than +, which means, in the above example, the `1+var1+var2`
 * term would construct a single, compound variable rather than adding 3 separate variables to the
 * equation.
 */
class Equation {
    public:
        /// Not default constructible
        Equation() = delete;

        /// Creates an equation with the given dependent variable
        Equation(const Variable &y);
        /// Creates an equation with the given dependent variable, accepting a rvalue reference
        template <class V, typename = typename std::enable_if<
            std::is_base_of<Variable, V>::value and std::is_move_constructible<V>::value>::type>
        Equation(V &&y);

        // Forward declaration
        class Proxy;

        /** Adds a variable (lvalue version).  If the value is a ConstantVariable with value 0, the
         * constant is removed from the equation; other ConstantVariable values replace (or re-add,
         * if previously removed) the constant.  Other types of variables are added to the equation.
         */
        Proxy operator % (const Variable &var);

        /** Adding a temporary variable works like above, but first stores the temporary by moving
         * it.  Only move constructible types are permitted.
         */
        template <class V, typename = typename std::enable_if<
            std::is_base_of<Variable, V>::value and std::is_move_constructible<V>::value>::type>
        Proxy operator % (V &&var);

        /// `% d` for double d is equivalent to `% ConstantVariable(d)`.
        Proxy operator % (double c);

        /// Accesses the dependent variable
        const Variable& depVar() const;

        /// Returns the number of independent variables.
        unsigned int numVars() const;

        /** Returns true if the model includes a constant, false if it does not.  If the model has a
         * constant, it is guaranteed to be the first Variable accessed via begin().
         */
        bool hasConstant() const;

        /** Returns the names of the variables in the same order they will be returned by
         * begin()/end().
         */
        std::vector<std::string> names() const;

        /** Const access to the beginning iterator of independent variables.  If an explicit
         * constant of 0 has been added to the model, this will skip the constant; otherwise the
         * constant will be the first element.
         */
        const std::list<std::reference_wrapper<const Variable>>::const_iterator begin() const;

        /** Const access to the past-the-end iterator of independent variables.
         */
        const std::list<std::reference_wrapper<const Variable>>::const_iterator end() const;

        /** Overloaded so that an Equation can be sent to an output stream, resulting in output such
         * as `y ~ const + x1 + x2`.
         */
        friend std::ostream& operator<<(std::ostream &os, const Equation &eq);

    private:
        // If given rvalue reference, move and store them here
        std::list<std::shared_ptr<Variable>> private_vars_;

        /// Helper function to initialize the private_vars_ list with a single value
        template <class V>
        static std::list<std::shared_ptr<Variable>> private_vars_initializer(V &&v);

        // The constant; if called with a new ConstantVariable, we replace this one.  If the
        // constant value becomes 0, independent variables list iteration won't include the
        // constant.  Defaults to 1.
        ConstantVariable const_;

    protected:
        /// The dependent variable
        const Variable &dep_var_;

        /// The independent variables; the first is always a ConstantVariable
        std::list<std::reference_wrapper<const Variable>> indep_vars_{{const_}};

        /// Internal method to add a variable to the model
        void addVar(const Variable &var);

        /** Internal method to add a variable to the model; this moves and stores the given
         * variable, then calls the lvalue version of the method.
         */
        template <class V, typename = typename std::enable_if<
            std::is_base_of<Variable, V>::value and std::is_move_constructible<V>::value>::type>
        void addVar(V &&var);
};

/// Proxy object returned by << that allows more variables to be appended with +
class Equation::Proxy final {
    public:
        /// Appends another variable to the model; see Equation::operator<<
        Proxy& operator + (const Variable &var);
        /// Appends another variable to the model; see Equation::operator<<
        template <class V, typename = typename std::enable_if<
            std::is_base_of<Variable, V>::value and std::is_move_constructible<V>::value>::type>
        Proxy& operator + (V &&var);
    private:
        Equation &eq_;
        Proxy(Equation &eq);
        friend class Equation;
};

template <class V>
std::list<std::shared_ptr<Variable>> Equation::private_vars_initializer(V &&v) {
    std::list<std::shared_ptr<Variable>> lst;
    lst.emplace_back(new V(std::move(v)));
    return lst;
}

template <class V, typename>
Equation::Equation(V &&y) : private_vars_(private_vars_initializer(std::move(y))), dep_var_(*private_vars_.front())
{}

template <class V, typename>
Equation::Proxy Equation::operator % (V &&var) {
    addVar(std::move(var));
    return Proxy(*this);
}

template <class V, typename>
Equation::Proxy& Equation::Proxy::operator + (V &&var) {
    eq_.addVar(std::move(var));
    return *this;
}

/// Specialization for a ConstantVariable
template <>
void Equation::addVar<ConstantVariable>(ConstantVariable &&v);

template <class V, typename>
void Equation::addVar(V &&v) {
    private_vars_.emplace_back(new V(std::move(v)));
    addVar(*private_vars_.back());
}

}}
