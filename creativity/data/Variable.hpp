#pragma once
#include <cmath>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <Eigen/Core>

namespace creativity { namespace data {

/// Helper wrapper around std::to_string that removes insignificant 0s and a .
std::string to_string(double d);

/** Abstract base class for model variables.  Subclasses of this are used for simple variables,
 * multiplications of variables, etc.
 */
class Variable {
    public:
        /// Virtual destructor
        virtual ~Variable() = default;

        /// Returns a name for this variable
        virtual std::string name() const = 0;

        /** Returns a name for this variable, but with surrounding brackets if this variable needs
         * it--for example, for a multiplication of two variables.  For simple variables and
         * constants, and some non-simple variables such as Logarithm, the brackets are not needed.
         *
         * \param bracketL the left bracket; defaults to a left parenthesis
         * \param bracketR the right bracket; defaults to a right parenthesis
         */
        virtual std::string nameBracketed(const std::string &bracketL = "(", const std::string &bracketR = ")") const;

        /** Populates the given column reference with this variable's values.
         *
         * The `offset` and `trim` parameters allow you to exclude leading or trailing elements from
         * the source data.  Note that these elements are not excluded from the populated column:
         * for example, when the data is backed by a column, `source.size() - offset - trim` must
         * equal `column.size()`.
         *
         * \param column the column to populate
         * \param offset if specified and non-zero, the number of values to skip from the
         * beginning of the source data (but NOT `column`).
         * \param trim if specified and non-zero, the number of values to leave off the end of the
         * source data (but NOT `column`).
         * \throws Variable::SizeError if the given column's size, after adjusting for `offset` and
         * `trim`, is not compatible with this variable's size.
         */
        virtual void populate(Eigen::Ref<Eigen::VectorXd> column, unsigned int offset = 0, unsigned int trim = 0) const = 0;

        /** Shortcut wrapper around populate that creates a new column of the given size, calls
         * populate() with it, then returns it.
         *
         * If rows equals 0 (the default if not specified), the size is determined by a call to
         * size().  Note, however, that size() will throw an exception if called on a Variable
         * without an implied size, typically a ConstantVariable.
         */
        Eigen::VectorXd values(unsigned int rows = 0, unsigned int offset = 0, unsigned int trim = 0) const;

        /** Returns the intrinsic size of this data.  If this Variable has no such intrinsic size,
         * such as a ConstantVariable, this throws a SizeError exception.  If this variable is an
         * expression with multiple sizes, this method return can return any of the sizes.
         */
        virtual unsigned int size() const = 0;

        /** Exception class thrown when attempting to populate a column from a source of an
         * incompatible or unknown size.
         */
        class SizeError : public std::logic_error {
            public:
                /// Constructs with default message about populate() size mismatch
                SizeError();
                /// Constructs with custom message
                SizeError(const std::string &why);
        };
};

/** Special Variable subclass for representing a constant value.  When populating, the constant will
 * be copied into every requested position of the given column.
 */
class ConstantVariable final : public Variable {
    public:
        /// Default constructor creates a constant with the value 1.0
        ConstantVariable() = default;
        /** Creates a constant with the given value (this constructor also allows implicit
         * conversion of a double to a ConstantVariable).
         */
        ConstantVariable(double c);
        /// Copies the constant into the given column positions.
        virtual void populate(Eigen::Ref<Eigen::VectorXd> column, unsigned int offset = 0, unsigned int trim = 0) const override;

        /** Returns a name.  If the constant is 1.0, returns "const"; if 0.0, returns "noconst";
         * otherwise returns `to_string(constant)`.
         */
        virtual std::string name() const override;

        /// Overridden to not bracket the name
        virtual std::string nameBracketed(const std::string& = "(", const std::string& = ")") const override;

        /// Accesses the constant value
        const double& value() const;

        /// Throws a SizeError exception (a ConstantVariable has no intrinsic size)
        unsigned int size() const override;

    private:
        /// The actual constant value
        double c_ = 1.0;
};

/** Wrapper class around a simple column, where values are exactly the value in the column. */
class SimpleVariable : public Variable {
    public:
        /// Not default constructible
        SimpleVariable() = delete;
        /// Constructs 
        SimpleVariable(const std::string &name, const Eigen::Ref<const Eigen::VectorXd> &values);
        /// Copies the source column values into `column`
        virtual void populate(Eigen::Ref<Eigen::VectorXd> column, unsigned int offset = 0, unsigned int trim = 0) const override;
        /// Returns the name given during construction
        virtual std::string name() const override;
        /// Overridden to not bracket the name
        virtual std::string nameBracketed(const std::string& = "(", const std::string& = ")") const override;
        /// Returns the size of the underlying data column
        unsigned int size() const override;

    protected:
        /// The name of this variable
        const std::string name_;
        /// A reference to the stored column of values
        Eigen::Ref<const Eigen::VectorXd> col_;
};

/** Returns true if V is a simple variable type, that is, a value that only requires copying but not
 * other operations.  This is primarily used to determine when parentheses need to be added to the
 * names of complex types.  The default implementation returns false; classes that are simple types
 * should provide specializations that return true.
 */
template <class V> constexpr bool is_simple() { return false; }
/// Specialization of is_simple for SimpleVariable that returns true
template <> constexpr bool is_simple<SimpleVariable>() { return true; }
/// Specialization of is_simple for ConstantVariable that returns true
template <> constexpr bool is_simple<ConstantVariable>() { return true; }

/// Common base class for all Addition<A,B> classes
class AdditionBase : public Variable {};
/// Common base class for all Multiplication<A,B> classes
class MultiplicationBase : public Variable {};
/// Common base class for all Division<A,B> classes
class DivisionBase : public Variable {};
/// Common base class for all Power<V> classes
class PowerBase : public Variable {};
/// Common base class for all Exponential<V> classes
class ExponentialBase : public Variable {};
/// Common base class for all Logarithm<V> classes
class LogarithmBase : public Variable {};

/** Wrapper class around a multiplication of the elements of two Variable subclassess.  This can be
 * constructed explicitly or implicitly by the multiplication operator of two Variable objects.
 */
template <class Left, class Right, typename = typename std::enable_if<std::is_base_of<Variable, Left>::value and std::is_base_of<Variable, Right>::value>::type>
class Multiplication : public MultiplicationBase {
    public:
        // Default-constructible iff Left and Right are default-constructible

        /** Multiplies two Variables together.  (Note that the multiplication is only done when
         * populate() is called to determine the actual value.)
         */
        Multiplication(Left left, Right right) : left_{std::move(left)}, right_{std::move(right)} {}

        /// Calculates the component-wise multiplication of the two variables.
        virtual void populate(Eigen::Ref<Eigen::VectorXd> column, unsigned int offset = 0, unsigned int trim = 0) const override {
            left_.populate(column, offset, trim);
            column.array() *= right_.values(column.size(), offset, trim).array();
        }

        /** Returns the name by joining together left and right names with "*".  If one or the other
         * is something other than a SimpleVariable or ConstantVariable, the name is surrounded with
         * parentheses.
         *
         * As a special case, if Left is a ConstantVariable equal to -1, and Right is not a
         * ConstantVariable, this results in "-name" instead of "-1*name".
         */
        virtual std::string name() const override {
            // Attempt some bracket collapses where order-of-operations allows
            std::string left(
                    std::is_base_of<MultiplicationBase, Left>::value or
                    std::is_base_of<DivisionBase, Left>::value or
                    std::is_base_of<PowerBase, Left>::value or
                    std::is_base_of<ExponentialBase, Left>::value or
                    std::is_base_of<LogarithmBase, Left>::value
                    ? left_.name() : left_.nameBracketed());
            std::string op("*");
            std::string right(
                    std::is_base_of<MultiplicationBase, Right>::value or
                    std::is_base_of<DivisionBase, Right>::value or
                    std::is_base_of<PowerBase, Right>::value or
                    std::is_base_of<ExponentialBase, Right>::value or
                    std::is_base_of<LogarithmBase, Right>::value
                    ? right_.name() : right_.nameBracketed());

            if (std::is_same<Left, ConstantVariable>::value and not std::is_same<Right, ConstantVariable>::value
                    // This dynamic_cast isn't doing anything since the first is_same above is true
                    and dynamic_cast<const ConstantVariable&>(left_).value() == -1.0) {
                left.clear(); op = "-";
            }
            else if (std::is_same<Right, ConstantVariable>::value and not std::is_same<Left, ConstantVariable>::value
                    and dynamic_cast<const ConstantVariable&>(right_).value() == -1.0) {
                right = left;
                op = "-";
                left.clear();
            }

            return left + op + right;
        }

        /// Examines the data set to attempt to find and return a data size.
        unsigned int size() const override {
            try { return left_.size(); }
            catch (SizeError&) {
                // If both sides throw, we're a complex operation of just constants (so we want to
                // throw):
                return right_.size();
            }
        }


    protected:
        /// The left-hand Variable
        Left left_;
        /// The right-hand Variable
        Right right_;
};

/** Multiplies two Variable subclass object together, returning a new Multiplication object. */
template <class Left, class Right, typename = typename std::enable_if<std::is_base_of<Variable, Left>::value and std::is_base_of<Variable, Right>::value>::type>
Multiplication<Left, Right> operator* (const Left &left, const Right &right) {
    return Multiplication<Left, Right>(left, right);
}

/// Multiplies a Variable by a constant.
template <class V, typename = typename std::enable_if<std::is_base_of<Variable, V>::value>::type>
Multiplication<ConstantVariable, V> operator* (const V &v, double c) {
    return Multiplication<ConstantVariable, V>(c, v);
}
/// Multiplies a constant by a Variable.
template <class V, typename = typename std::enable_if<std::is_base_of<Variable, V>::value>::type>
Multiplication<ConstantVariable, V> operator* (double c, const V &v) {
    return Multiplication<ConstantVariable, V>(c, v);
}
/// Unary negation of a Variable is converted to a multiplication by -1
template <class V, typename = typename std::enable_if<std::is_base_of<Variable, V>::value>::type>
Multiplication<ConstantVariable, V> operator- (const V &v) {
    return -1.0 * v;
}

/** Wrapper class around a addition of the elements of two Variable subclassess.  This can be
 * constructed explicitly or implicitly by the addition or subtraction operator of two Variable
 * objects.
 */
template <class Left, class Right, typename = typename std::enable_if<std::is_base_of<Variable, Left>::value and std::is_base_of<Variable, Right>::value>::type>
class Addition : public AdditionBase {
    public:
        // Default-constructible iff Left and Right are default-constructible

        /** Adds two Variables together.  (Note that the addition is only done when populate() is
         * called to determine the actual value.)
         */
        Addition(Left left, Right right) : left_{std::move(left)}, right_{std::move(right)} {}

        /// Calculates the component-wise addition of the two variables.
        virtual void populate(Eigen::Ref<Eigen::VectorXd> column, unsigned int offset = 0, unsigned int trim = 0) const override {
            left_.populate(column, offset, trim);
            column.array() += right_.values(column.size(), offset, trim).array();
        }

        /** Returns the name by joining together left and right names with "+".  Since there is
         * nothing lower in the order or operations, this returns unbracketed left- and right-hand
         * side operands.
         */
        virtual std::string name() const override { return left_.name() + "+" + right_.name(); }

        /// Examines the data set to attempt to find and return a data size.
        unsigned int size() const override {
            try { return left_.size(); }
            catch (SizeError&) {
                // If both sides throw, we're a complex operation of just constants (so we want to
                // throw):
                return right_.size();
            }
        }

    protected:
        /// The left-hand Variable
        Left left_;
        /// The right-hand Variable
        Right right_;
};

/** Adds two Variable subclass object together, returning a new Addition object. */
template <class Left, class Right, typename = typename std::enable_if<std::is_base_of<Variable, Left>::value and std::is_base_of<Variable, Right>::value>::type>
Addition<Left, Right> operator+ (const Left &left, const Right &right) {
    return Addition<Left, Right>(left, right);
}

/// Adds a constant to a Variable; the constant is converted to a ConstantVariable
template <class V, typename = typename std::enable_if<std::is_base_of<Variable, V>::value>::type>
Addition<V, ConstantVariable> operator+ (const V &v, double c) {
    return Addition<V, ConstantVariable>(v, c);
}
/// Adds a Variable to a constant; the constant is converted to a ConstantVariable
template <class V, typename = typename std::enable_if<std::is_base_of<Variable, V>::value>::type>
Addition<ConstantVariable, V> operator+ (double c, const V &v) {
    return Addition<ConstantVariable, V>(c, v);
}
/// Subtracts one Variable from another; the subtracted value is multiplied by -1.
template <class Left, class Right, typename = typename std::enable_if<std::is_base_of<Variable, Left>::value and std::is_base_of<Variable, Right>::value>::type>
Addition<Left, Multiplication<ConstantVariable, Right>> operator- (const Left &left, const Right &right) {
    return Addition<Left, Multiplication<ConstantVariable, Right>>(left, -right);
}
/// Subtracts a constant from a Variable; the constant is converted to a ConstantVariable
template <class V, typename = typename std::enable_if<std::is_base_of<Variable, V>::value>::type>
Addition<V, ConstantVariable> operator- (const V &v, double c) {
    return Addition<V, ConstantVariable>(v, -c);
}
/** Subtracts a Variable from a constant; the constant is converted to a ConstantVariable, and the
 * Variable is multiplied by -1.
 */
template <class V, typename = typename std::enable_if<std::is_base_of<Variable, V>::value>::type>
Addition<ConstantVariable, Multiplication<ConstantVariable, V>> operator- (double c, const V &v) {
    return Addition<ConstantVariable, Multiplication<ConstantVariable, V>>(c, -v);
}

/** Class that returns the coefficient-wise division of elements in one Variable by corresponding
 * elements in the other Variable.
 *
 * The class can be used directly or via the overload of the '/' operator.
 */
template <class Top, class Bottom, typename = typename std::enable_if<std::is_base_of<Variable, Top>::value and std::is_base_of<Variable, Bottom>::value>::type>
class Division : public DivisionBase {
    public:
        // Default-constructible iff Top and Bottom are default-constructible

        /** Divides one Variable by another.  (Note that the division is only done when populate()
         * is called to determine the actual value.)
         */
        Division(Top numerator, Bottom denominator) : top_{std::move(numerator)}, bottom_{std::move(denominator)} {}

        /// Calculates the component-wise multiplication of the two variables.
        virtual void populate(Eigen::Ref<Eigen::VectorXd> column, unsigned int offset = 0, unsigned int trim = 0) const override {
            top_.populate(column, offset, trim);
            column.array() /= bottom_.values(column.size(), offset, trim).array();
        }

        /** Returns the name by joining together left and right names with "/".  If one or the other
         * is something other than a SimpleVariable or ConstantVariable, the name is surrounded with
         * parentheses as needed.
         */
        virtual std::string name() const override {
            // This is similar (but not identical) to the Multiplication version: the main
            // difference is in the denominator, which has different precedence rules than the top.

            // Attempt some bracket collapses where order-of-operations allows
            std::string left(
                    std::is_base_of<MultiplicationBase, Top>::value or
                    std::is_base_of<DivisionBase, Top>::value or
                    std::is_base_of<PowerBase, Top>::value or
                    std::is_base_of<ExponentialBase, Top>::value or
                    std::is_base_of<LogarithmBase, Top>::value
                    ? top_.name() : top_.nameBracketed());
            std::string op("/");
            std::string right(
                    std::is_base_of<PowerBase, Bottom>::value or
                    std::is_base_of<ExponentialBase, Bottom>::value or
                    std::is_base_of<LogarithmBase, Bottom>::value
                    ? bottom_.name() : bottom_.nameBracketed());

            if (std::is_same<Top, ConstantVariable>::value and not std::is_same<Bottom, ConstantVariable>::value
                    // This dynamic_cast isn't doing anything since the first is_same above is true
                    and dynamic_cast<const ConstantVariable&>(top_).value() == -1.0) {
                left.clear(); op = "-";
            }
            else if (std::is_same<Bottom, ConstantVariable>::value and not std::is_same<Top, ConstantVariable>::value
                    and dynamic_cast<const ConstantVariable&>(bottom_).value() == -1.0) {
                right = left;
                op = "-";
                left.clear();
            }

            return left + op + right;
        }

        /// Examines the data set to attempt to find and return a data size.
        unsigned int size() const override {
            try { return top_.size(); }
            catch (SizeError&) {
                // If both sides throw, we're a complex operation of just constants (so we want to
                // throw):
                return bottom_.size();
            }
        }

    protected:
        /// The left-hand Variable
        Top top_;
        /// The right-hand Variable
        Bottom bottom_;
};

/** Divides one Variable by another, returning a new Division object. */
template <class Top, class Bottom, typename = typename std::enable_if<std::is_base_of<Variable, Top>::value and std::is_base_of<Variable, Bottom>::value>::type>
Division<Top, Bottom> operator/ (const Top &numerator, const Bottom &denominator) {
    return Division<Top, Bottom>(numerator, denominator);
}
/// Divides a Variable by a constant.
template <class V, typename = typename std::enable_if<std::is_base_of<Variable, V>::value>::type>
Division<V, ConstantVariable> operator/ (const V &v, double c) {
    return Division<V, ConstantVariable>(v, c);
}
/// Divides a constant by a Variable.
template <class V, typename = typename std::enable_if<std::is_base_of<Variable, V>::value>::type>
Division<ConstantVariable, V> operator/ (double c, const V &v) {
    return Division<ConstantVariable, V>(c, v);
}

/** Raises a Variable's values to a power.  This is done using std::pow unless the power is one of
 * the special values -1, 0.5, 1, 2, or 3.
 *
 * The class can be used explicitly, but a specialization of `std::pow(var, power)` and an overload
 * of `var ^ double` are also available.
 */
template <class V, typename = typename std::enable_if<std::is_base_of<Variable, V>::value>::type>
class Power : public PowerBase {
    public:
        /// Not default constructible
        Power() = delete;

        /** Wraps around a Variable to provide exponentiation of the variable's values. */
        Power(V var, double power) : var_{std::move(var)}, power_{power} {}

        /// Calculates and stores the raised values
        virtual void populate(Eigen::Ref<Eigen::VectorXd> column, unsigned int offset = 0, unsigned int trim = 0) const override {
            var_.populate(column, offset, trim);
            if (power_ == 2)
                column = column.array().square();
            else if (power_ == 3)
                column = column.array().cube();
            else if (power_ == 0.5)
                column = column.array().sqrt();
            else if (power_ == -1)
                column = column.array().inverse();
            else if (power_ != 1)
                column = column.array().pow(power_);
        }

        /** Returns the concatenation the underlying Variable name() with "^" and the power.  If V
         * is something other than SimpleVariable or ConstantVariable, the underlying name is also
         * surrounded by parentheses.
         */
        virtual std::string name() const override { return var_.nameBracketed() + "^" + to_string(power_); }

        /// Examines the data set to attempt to find and return a data size.
        unsigned int size() const override { return var_.size(); }

    protected:
        /// The Variable subclass with values to be inverted
        V var_;
        /// The power to which to raise the variable
        double power_;
};

/** `Variable ^ power` returns a Power object. */
template <class V, typename = typename std::enable_if<std::is_base_of<Variable, V>::value>::type>
Power<V> operator^ (const V &val, double pow) {
    return Power<V>(val, pow);
}

/** Raises a base value to a Variable's value.  This uses std::pow for exponentiation except in the
 * special cases where `base == std::exp(1)` and `base == 2`, where std::exp and std::exp2 are used
 * instead.
 *
 * The class can be used explicitly, but a specializations of `std::pow(base, var)`,
 * `std::exp(var)`, and `std::exp2(var)` and an overload of `double ^ var` are also available.
 */
template <class V, typename = typename std::enable_if<std::is_base_of<Variable, V>::value>::type>
class Exponential : public ExponentialBase {
    public:
        /// Not default constructible
        Exponential() = delete;

        /** Wraps around a Variable to provide exponentiation of the variable's values. */
        Exponential(double base, V var) : var_{std::move(var)}, base_{base} {}

        /** Wraps around a Variable to provide exponentiation of the variable's values using Euler's
         * number (e).
         */
        Exponential(V var) : var_{std::move(var)}, base_{std::exp(1)} {}

        /// Calculates and stores the exponential values
        virtual void populate(Eigen::Ref<Eigen::VectorXd> column, unsigned int offset = 0, unsigned int trim = 0) const override {
            var_.populate(column, offset, trim);
            if (base_ == std::exp(1))
                column = column.array().exp().matrix();
            else if (base_ == 2)
                column = column.unaryExpr([](double c) { return std::exp2(c); });
            else
                column = column.unaryExpr([this](double c) { return std::pow(base_, c); });
        }

        /** Returns a string representation of the variable.  If the base equals `std::exp(1)`, this
         * is "exp(name)"; otherwise the name is "base^name", with parentheses added to name if it
         * is a complex variable.
         */
        virtual std::string name() const override {
            if (base_ == std::exp(1)) return "exp(" + var_.name() + ")"; // NB: don't need bracketed name here
            else return to_string(base_) + "^" + var_.nameBracketed();
        }

        /// Returns the name of the variable such as '{2^name}' or 'exp(name)' or '{2^log(name)}'
        virtual std::string nameBracketed(const std::string &bracketL = "(", const std::string &bracketR = ")") const override {
            if (base_ == std::exp(1)) return "exp(" + var_.name() + ")";
            else return bracketL + to_string(base_) + "^" + var_.nameBracketed() + bracketR;
        }

        /// Examines the data set to attempt to find and return a data size.
        unsigned int size() const override { return var_.size(); }

    protected:
        /// The Variable subclass with values to be inverted
        V var_;
        /// The power to which to raise the variable
        double base_;
};

/** `base ^ Variable` returns an Exponential object. */
template <class V, typename = typename std::enable_if<std::is_base_of<Variable, V>::value>::type>
Exponential<V> operator^ (double base, const V &val) {
    return Exponential<V>(base, val);
}

/** Takes a natural logarithm of a Variable's values.
 */
template <class V, typename = typename std::enable_if<std::is_base_of<Variable, V>::value>::type>
class Logarithm : public LogarithmBase {
    public:
        /// Not default constructible
        Logarithm() = delete;

        /** Wraps around a Variable to provide a natural logarithm calculation of the variable's
         * values.
         */
        Logarithm(V var) : var_{std::move(var)} {}

        /// Calculates and stores the logarithm values
        virtual void populate(Eigen::Ref<Eigen::VectorXd> column, unsigned int offset = 0, unsigned int trim = 0) const override {
            var_.populate(column, offset, trim);
            column = column.array().log().matrix();
        }

        /** Returns a string representation of the variable, which is `log(name)`, where name is the
         * name of the variable being log'ed.
         */
        virtual std::string name() const override {
            return "log(" + var_.name() + ")";
        }

        virtual std::string nameBracketed(const std::string& = "(", const std::string& = ")") const override {
            return name();
        }

        /// Examines the data set to attempt to find and return a data size.
        unsigned int size() const override { return var_.size(); }

    protected:
        /// The Variable subclass with values to be log'ed
        V var_;
};

}}

namespace std {
/// std::exp specialization for a Variable.  Returns an Exponential variable wrapper.
template <class V>
typename std::enable_if<std::is_base_of<creativity::data::Variable, V>::value, creativity::data::Exponential<V>>::type
exp(const V &var) { return var; }
/// std::exp2 specialization for a Variable.  Returns an Exponential variable wrapper with base 2.
template <class V>
typename std::enable_if<std::is_base_of<creativity::data::Variable, V>::value, creativity::data::Exponential<V>>::type
exp2(const V &var) { return creativity::data::Exponential<V>(2, var); }
/// std::log specialization for a Variable.  Returns a Logarithm variable wrapper.
template <class V>
typename std::enable_if<std::is_base_of<creativity::data::Variable, V>::value, creativity::data::Logarithm<V>>::type
log(const V &var) { return var; }
/// std::sqrt specialization for a Variable.  Returns a Power variable wrapper with power = 0.5.
template <class V>
typename std::enable_if<std::is_base_of<creativity::data::Variable, V>::value, creativity::data::Power<V>>::type
sqrt(const V &var) { return creativity::data::Power<V>(var, 0.5); }
/// std::pow specialization for a Variable raised to a numeric power.  Returns a Power variable wrapper.
template <class V>
typename std::enable_if<std::is_base_of<creativity::data::Variable, V>::value, creativity::data::Power<V>>::type
pow(const V &var, double power) { return creativity::data::Power<V>(var, power); }
/// std::pow specialization for a double raised to a Variable.  Returns an Exponential variable wrapper.
template <class V>
typename std::enable_if<std::is_base_of<creativity::data::Variable, V>::value, creativity::data::Exponential<V>>::type
pow(double base, const V &var) { return creativity::data::Exponential<V>(base, var); }
}
