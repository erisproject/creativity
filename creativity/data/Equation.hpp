#pragma once
#include <list>
#include <memory>
#include <functional>
#include <utility>
#include <ostream>
#include <string>
#include <type_traits>
#include <vector>
#include "creativity/data/Variable.hpp"

namespace creativity { namespace data {

/** Class to store a equation.  The class is constructed with the dependent variable, then has an
 * equation added to it using the % operator.  For example:
 *
 *     // Constructing a single equation
 *     Equation model(depvar);
 *     model % 1 + var1 + var2;
 *     model % var3;
 *
 *     // Building on another equation:
 *     Equation model2 = model1 + var4;
 *
 *     // rvalue construction
 *     Equation(depvar) + 1 + var1 + var2;
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

        /// Move constructor
        Equation(Equation &&move) = default;

        /// Copy constructor
        Equation(const Equation &copy) = default;

        /// Creates an equation with the given dependent variable
        Equation(const std::shared_ptr<const Variable> &y);

        // Forward declaration
        class Proxy;

        /** Adds a variable (lvalue version).  If the value is a ConstantVariable with value 0, the
         * constant is removed from the equation; other ConstantVariable values replace (or re-add,
         * if previously removed) the constant.  Other types of variables are added to the equation.
         */
        Proxy operator % (const std::shared_ptr<const Variable> &var);

        /// `% d` for double d is equivalent to `% ConstantVariable(d)`.
        Proxy operator % (double c);

        /** The addition operator duplicates the called-upon Equation an adds a new term to the
         * duplicate, then returns it.
         */
        Equation operator + (const std::shared_ptr<const Variable> &var) const &;

        /** The addition operator called on a temporary adds a new term to the temporary, then
         * returns it.
         */
        Equation operator + (const std::shared_ptr<const Variable> &var) &&;

        /// Adding a double constant (typically 0 or 1) converts the constant to a ConstantVariable
        Equation operator + (double c) const &;

        /// Adding a double constant (typically 0 or 1) converts the constant to a ConstantVariable
        Equation operator + (double c) &&;

        /// Accesses the dependent variable
        std::shared_ptr<const Variable> depVar() const;

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
        std::list<std::shared_ptr<const Variable>>::const_iterator begin() const;

        /** Const access to the past-the-end iterator of independent variables.
         */
        std::list<std::shared_ptr<const Variable>>::const_iterator end() const;

        /** Overloaded so that an Equation can be sent to an output stream, resulting in output such
         * as `y ~ const + x1 + x2`.
         */
        friend std::ostream& operator<<(std::ostream &os, const Equation &eq);

    protected:
        /// The dependent variable
        std::shared_ptr<const Variable> dep_var_;

        /// The independent variables; the first is always a ConstantVariable
        std::list<std::shared_ptr<const Variable>> indep_vars_{{ConstantVariable::create()}};

        /// Internal method to add a variable to the model
        void addVar(const std::shared_ptr<const Variable> &var);
};

/// Proxy object returned by << that allows more variables to be appended with +
class Equation::Proxy final {
    public:
        /// Appends another variable to the model; see Equation::operator%
        Proxy& operator + (const std::shared_ptr<const Variable> &var);
    private:
        Equation &eq_;
        Proxy(Equation &eq);
        friend class Equation;
};

}}
