#pragma once
#include <cmath>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <memory>
#include <Eigen/Core>

namespace creativity { namespace data {

/// Helper wrapper around std::to_string that removes insignificant 0s and a .
std::string to_string(double d);

#define CREATE_SHARED_WRAPPER(C) template <class... Args> static std::shared_ptr<C> create(Args... args) { return std::shared_ptr<C>(new C(std::forward<Args>(args)...)); }

/** Abstract base class for model variables.  Subclasses of this are used for simple variables,
 * multiplications of variables, etc.
 */
class Variable : public std::enable_shared_from_this<Variable> {
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
    protected:
        /// Default constructor creates a constant with the value 1.0
        ConstantVariable() = default;
        /** Creates a constant with the given value (this constructor also allows implicit
         * conversion of a double to a ConstantVariable).
         */
        ConstantVariable(double c);
    public:
        /// Forwards arguments to the protected constructor and returns a shared_ptr to the created object.
        CREATE_SHARED_WRAPPER(ConstantVariable)

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
    protected:
        /// Not default constructible
        SimpleVariable() = delete;

        /** Constructs a SimpleVariable that copies values from the given Vector.
         */
        SimpleVariable(const std::string &name, const Eigen::Ref<const Eigen::VectorXd> values);

    public:
        /// Forwards arguments to the protected constructor and returns a shared_ptr to the created object.
        CREATE_SHARED_WRAPPER(SimpleVariable)

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
        /// The column of values
        Eigen::VectorXd col_;
};

/** Base class for composite variables with two Variable components. */
class BinaryExpr : public Variable {
    public:
        /** Returns the size of the binary expression, if either left or right variables have a
         * size.  If neither do, throws a SizeError exception.
         */
        virtual unsigned int size() const override;

    protected:
        /// Not default-constructible
        BinaryExpr() = delete;

        /// Constructs a binary variable with the left and right variables.
        BinaryExpr(const std::shared_ptr<const Variable> &left, const std::shared_ptr<const Variable> &right);

        /// The left-hand Variable
        std::shared_ptr<const Variable> left_;
        /// The right-hand Variable
        std::shared_ptr<const Variable> right_;
};

/** Base class for composite variables with a single variable.  This isn't always strictly a unary
 * expression--it can also be used by classes that, for example, also handle a double value. */
class UnaryExpr : public Variable {
    public:
        /** Returns the size of the unary expression, if the variable has a size.  If it doesn't
         * throws a SizeError exception.
         */
        virtual unsigned int size() const override;

    protected:
        /// Not default-constructible
        UnaryExpr() = delete;

        /// Constructs a unary variable
        UnaryExpr(const std::shared_ptr<const Variable> &var);

        /// The underlying unary variable
        std::shared_ptr<const Variable> var_;
};


/** Wrapper class around a multiplication of the elements of two Variables.  This can be constructed
 * explicitly or implicitly by the multiplication operator of two Variable objects.
 */
class Multiplication : public BinaryExpr {
    protected:
        /// Inherit constructor(s) from BinaryExpr
        using BinaryExpr::BinaryExpr;

    public:
        /// Forwards arguments to the protected constructor and returns a shared_ptr to the created object.
        CREATE_SHARED_WRAPPER(Multiplication)

        /// Calculates the component-wise multiplication of the two variables.
        virtual void populate(Eigen::Ref<Eigen::VectorXd> column, unsigned int offset = 0, unsigned int trim = 0) const override;

        /** Returns the name by joining together left and right names with "*".  If one or the other
         * is something other than a SimpleVariable or ConstantVariable, the name is surrounded with
         * parentheses.
         *
         * As a special case, if the left variable is a ConstantVariable equal to -1, and Right is not a
         * ConstantVariable, this results in "-name" instead of "-1*name".
         */
        virtual std::string name() const override;
};

/** Wrapper class around a addition of the elements of two Variables.  This can be constructed
 * explicitly or implicitly by the addition or subtraction operator of two Variable objects.
 */
class Addition : public BinaryExpr {
    protected:
        /// Inherit constructor(s) from BinaryExpr
        using BinaryExpr::BinaryExpr;

    public:
        /// Forwards arguments to the protected constructor and returns a shared_ptr to the created object.
        CREATE_SHARED_WRAPPER(Addition)

        /// Calculates the component-wise addition of the two variables.
        virtual void populate(Eigen::Ref<Eigen::VectorXd> column, unsigned int offset = 0, unsigned int trim = 0) const override;

        /** Returns the name by joining together left and right names with "+".  Since there is
         * nothing lower in the order or operations, this returns unbracketed left- and right-hand
         * side operands.
         */
        virtual std::string name() const override;
};

/** Class that returns the coefficient-wise division of elements in one Variable by corresponding
 * elements in the other Variable.
 *
 * The class can be used directly or via the overload of the '/' operator.
 */
class Division : public BinaryExpr {
    protected:
        /// Inherit constructor(s) from BinaryExpr
        using BinaryExpr::BinaryExpr;

    public:
        /// Forwards arguments to the protected constructor and returns a shared_ptr to the created object.
        CREATE_SHARED_WRAPPER(Division)

        /// Calculates the component-wise multiplication of the two variables.
        virtual void populate(Eigen::Ref<Eigen::VectorXd> column, unsigned int offset = 0, unsigned int trim = 0) const override;

        /** Returns the name by joining together left and right names with "/".  If one or the other
         * is something other than a SimpleVariable or ConstantVariable, the name is surrounded with
         * parentheses as needed.
         */
        virtual std::string name() const override;
};

/** Raises a Variable's values to a power.  This is done using std::pow unless the power is one of
 * the special values -1, 0.5, 1, 2, or 3.
 *
 * The class can be used explicitly, but a specialization of `std::pow(var, power)` and an overload
 * of `var ^ double` are also available.
 */
class Power : public UnaryExpr {
    protected:
        /// Not default constructible
        Power() = delete;

        /** Wraps around a Variable to provide exponentiation of the variable's values. */
        Power(const std::shared_ptr<const Variable> &var, double power);

    public:
        /// Forwards arguments to the protected constructor and returns a shared_ptr to the created object.
        CREATE_SHARED_WRAPPER(Power)

        /// Calculates and stores the raised values
        virtual void populate(Eigen::Ref<Eigen::VectorXd> column, unsigned int offset = 0, unsigned int trim = 0) const override;

        /** Returns the concatenation the underlying Variable name() with "^" and the power.  If V
         * is something other than SimpleVariable or ConstantVariable, the underlying name is also
         * surrounded by parentheses.
         */
        virtual std::string name() const override;

    protected:
        /// The power to which to raise the variable
        double power_;
};

/** Raises a base value to a Variable's value.  This uses std::pow for exponentiation except in the
 * special cases where `base == std::exp(1)` and `base == 2`, where std::exp and std::exp2 are used
 * instead.
 *
 * The class can be used explicitly, but a specializations of `std::pow(base, var)`,
 * `std::exp(var)`, and `std::exp2(var)` and an overload of `double ^ var` are also available.
 */
class Exponential : public UnaryExpr {
    protected:
        /** Wraps around a Variable to provide exponentiation of the variable's values. */
        Exponential(double base, const std::shared_ptr<const Variable> &var);

        /** Wraps around a Variable to provide exponentiation of the variable's values using Euler's
         * number (e).
         */
        Exponential(const std::shared_ptr<const Variable> &var);

    public:
        /// Forwards arguments to the protected constructor and returns a shared_ptr to the created object.
        CREATE_SHARED_WRAPPER(Exponential)

        /// Calculates and stores the exponential values
        virtual void populate(Eigen::Ref<Eigen::VectorXd> column, unsigned int offset = 0, unsigned int trim = 0) const override;

        /** Returns a string representation of the variable.  If the base equals `std::exp(1)`, this
         * is "exp(name)"; otherwise the name is "base^name", with parentheses added to name if it
         * is a complex variable.
         */
        virtual std::string name() const override;

        /// Returns the name of the variable such as '{2^name}' or 'exp(name)' or '{2^log(name)}'
        virtual std::string nameBracketed(const std::string &bracketL = "(", const std::string &bracketR = ")") const override;

    protected:
        /// The power to which to raise the variable
        double base_;
};

/** Takes a natural logarithm of a Variable's values.
 */
class Logarithm : public UnaryExpr {
    protected:
        /// Inherit constructors from UnaryExpr
        using UnaryExpr::UnaryExpr;

    public:
        /// Forwards arguments to the protected constructor and returns a shared_ptr to the created object.
        CREATE_SHARED_WRAPPER(Logarithm)

        /// Calculates and stores the logarithm values
        virtual void populate(Eigen::Ref<Eigen::VectorXd> column, unsigned int offset = 0, unsigned int trim = 0) const override;

        /** Returns a string representation of the variable, which is `log(name)`, where name is the
         * name of the variable being log'ed.
         */
        virtual std::string name() const override;

        virtual std::string nameBracketed(const std::string& = "(", const std::string& = ")") const override;
};

#undef CREATE_SHARED_WRAPPER

/** Multiplies two Variable objects together, returning a new Multiplication object. */
std::shared_ptr<Multiplication> operator* (const std::shared_ptr<const Variable> &left, const std::shared_ptr<const Variable> &right);
/// Multiplies a Variable by a constant.
std::shared_ptr<Multiplication> operator* (const std::shared_ptr<const Variable> &v, double c);
/// Multiplies a constant by a Variable.
std::shared_ptr<Multiplication> operator* (double c, const std::shared_ptr<const Variable> &v);
/// Unary negation of a Variable is converted to a multiplication by -1
std::shared_ptr<Multiplication> operator- (const std::shared_ptr<const Variable> &v);
/** Adds two Variable objects together, returning a new Addition object. */
std::shared_ptr<Addition> operator+ (const std::shared_ptr<const Variable> &left, const std::shared_ptr<const Variable> &right);
/// Adds a constant to a Variable; the constant is converted to a ConstantVariable
std::shared_ptr<Addition> operator+ (const std::shared_ptr<const Variable> &v, double c);
/// Adds a Variable to a constant; the constant is converted to a ConstantVariable
std::shared_ptr<Addition> operator+ (double c, const std::shared_ptr<const Variable> &v);
/// Subtracts one Variable from another; the subtracted value is multiplied by -1.
std::shared_ptr<Addition> operator- (const std::shared_ptr<const Variable> &left, const std::shared_ptr<const Variable> &right);
/// Subtracts a constant from a Variable; the constant is converted to a ConstantVariable
std::shared_ptr<Addition> operator- (const std::shared_ptr<const Variable> &v, double c);
/** Subtracts a Variable from a constant; the constant is converted to a ConstantVariable, and the
 * Variable is multiplied by -1.
 */
std::shared_ptr<Addition> operator- (double c, const std::shared_ptr<const Variable> &v);
/** Divides one Variable by another, returning a new Division object. */
std::shared_ptr<Division> operator/ (const std::shared_ptr<const Variable> &numerator, const std::shared_ptr<const Variable> &denominator);
/// Divides a Variable by a constant.
std::shared_ptr<Division> operator/ (const std::shared_ptr<const Variable> &v, double c);
/// Divides a constant by a Variable.
std::shared_ptr<Division> operator/ (double c, const std::shared_ptr<const Variable> &v);
/** `Variable ^ power` returns a Power object. */
std::shared_ptr<Power> operator^ (const std::shared_ptr<const Variable> &val, double pow);
/** `base ^ Variable` returns an Exponential object. */
std::shared_ptr<Exponential> operator^ (double base, const std::shared_ptr<const Variable> &val);

}}

namespace std {
/// std::exp specialization for a Variable.  Returns an Exponential variable wrapper.
shared_ptr<creativity::data::Exponential>
exp(const shared_ptr<const creativity::data::Variable> &var);
/// std::exp2 specialization for a Variable.  Returns an Exponential variable wrapper with base 2.
shared_ptr<creativity::data::Exponential>
exp2(const shared_ptr<const creativity::data::Variable> &var);
/// std::log specialization for a Variable.  Returns a Logarithm variable wrapper.
shared_ptr<creativity::data::Logarithm>
log(const shared_ptr<const creativity::data::Variable> &var);
/// std::sqrt specialization for a Variable.  Returns a Power variable wrapper with power = 0.5.
shared_ptr<creativity::data::Power>
sqrt(const shared_ptr<const creativity::data::Variable> &var);
/// std::pow specialization for a Variable raised to a numeric power.  Returns a Power variable wrapper.
shared_ptr<creativity::data::Power>
pow(const shared_ptr<const creativity::data::Variable> &var, double power);
/// std::pow specialization for a double raised to a Variable.  Returns an Exponential variable wrapper.
shared_ptr<creativity::data::Exponential>
pow(double base, const shared_ptr<const creativity::data::Variable> &var);
}
