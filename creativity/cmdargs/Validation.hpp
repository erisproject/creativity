#pragma once
#include "creativity/cmdargs/strings.hpp"
#include <boost/program_options/errors.hpp>
#include <string>

namespace creativity { namespace cmdargs {

/** Validation tag; any Validation class must (ultimately) inherit from this class.  This class does
 * nothing.
 */
class ValidationTag {};

/** Validation wrapper base class.  This base class does no actual value validation (except one, see
 * below).
 *
 * This should be inherited from virtually so that subclasses can inherit from multiple Validation
 * subclasses to enforce multiple validations at once.
 *
 * Note that there is one validation that occurs when using this: if `T` is an unsigned type,
 * validation of a Validation (or derived class) makes sure the given value is not negative; without
 * this, boost will cast negative input values to unsigned types, thus ending up with -1 becoming
 * (for an unsigned int) 4294967295.
 */
template <typename T>
class Validation : public ValidationTag {
    public:
        /// Constructs with an initial value
        Validation(T v) : val_(v) {}
        /// Implicit conversion to the stored value
        operator const T& () const { return val_; }
        /// The type T that this object validates
        using value_type = T;
        /// Returns a string representation of this validation object
        static std::string validationString() {
            if (std::is_unsigned<T>::value) { return type_string<T>() + u8"⩾0"; }
            return type_string<T>();
        }
        /// Virtual destructor
        virtual ~Validation() = default;

    protected:
        /// The stored value
        T val_;
};

/** Validation wrapper for options that have a minimum value.
 *
 * For technical reasons, the minimum must be given as a fraction of longs (non-type
 * template parameters cannot be non-integral types).
 *
 * \param T any numeric type.
 * \param min the minimum accepted value; if the value is fractional, this is the
 * numerator.
 * \param denom the denominator of the minimum value; defaults to 1.
 */
template <typename T, long min, long denom = 1>
class Min : public virtual Validation<T> {
    public:
        /// Constructor.  Throws if `v < min`.
        Min(T v) : Validation<T>(v) {
            if (*this < min / (double) denom)
                throw boost::program_options::validation_error(boost::program_options::validation_error::invalid_option_value);
        }

        /// Returns string representation of this validation
        static std::string validationString() { return type_string<T>() + u8"⩾" + (denom == 1 ? output_string(min) : output_string(min / (double) denom)); }
};

/** Validation wrapper for options that have a maximum value.
 *
 * For technical reasons, the minimum must be given as a fraction of longs (non-type
 * template parameters cannot be non-integral types).
 *
 * \param T any numeric type.
 * \param max the maximum accepted value; if fractional, this is the numerator.
 * \param denom the denominator of the maximum value; defaults to 1.
 */
template<typename T, long max, long denom = 1>
class Max : public virtual Validation<T> {
    public:
        /// Constructor.  Throws if `v > max`.
        Max(T v) : Validation<T>(v) {
            if (*this > max / (double) denom)
                throw boost::program_options::validation_error(boost::program_options::validation_error::invalid_option_value);
        }

        /// Returns string representation of this validation
        static std::string validationString() { return type_string<T>() + u8"⩽" + (denom == 1 ? output_string(max) : output_string(max / (double) denom)); }
};

/** Validation wrapper for options that have both a minimum and maximum value.
 *
 * \param T any numeric type
 * \param min the minimum accepted value; if fractional, this is the numerator.
 * \param max the maximum accepted value; if fractional, this is the numerator.
 * \param denom the denominator of both the min and max values; defaults to 1.
 */
template<typename T, long min, long max, long denom = 1>
class Range : public Min<T, min>, public Max<T, max> {
    public:
        /// Constructor.  Throws if `v < min` or `v > max`.
        Range(T v) : Validation<T>(v), Min<T, min, denom>(v), Max<T, max, denom>(v) {}

        /// Returns string representation of this validation
        static std::string validationString() { return (denom == 1 ? output_string(min) : output_string(min / (double) denom)) + u8"⩽" + type_string<T>() + u8"⩽" + (denom == 1 ? output_string(max) : output_string(max / (double) denom)); }
};

/** Validation wrapper for options that have a strict inequality boundary that the value
 * must be above, such as \f$v > 3\f$.
 *
 * \param T any numeric type
 * \param lower the lower bound that the value must be above (or the numerator of the value)
 * \param the denominator under `lower`.  Defaults to 1.
 */
template <typename T, long lower, long denom = 1>
class Above : public virtual Validation<T> {
    public:
        /// Constructor.  Throws if `v <= lower`.
        Above(T v) : Validation<T>(v) {
            if (*this <= lower / (double) denom)
                throw boost::program_options::validation_error(boost::program_options::validation_error::invalid_option_value);
        }

        /// Returns string representation of this validation
        static std::string validationString() { return type_string<T>() + u8">" + (denom == 1 ? output_string(lower) : output_string(lower / (double) denom)); }
};

/** Validation wrapper for options that have a strict inequality boundary that the value
 * must be below, such as \f$v < 3\f$.
 *
 * \param T any numeric type
 * \param upper the upper bound that the value must be above (or the numerator of the value)
 * \param the denominator under `lower`.  Defaults to 1.
 */
template <typename T, long upper, long denom = 1>
class Below : public virtual Validation<T> {
    public:
        /// Constructor.  Throws if `v >= upper`.
        Below(T v) : Validation<T>(v) {
            if (*this >= upper / (double) denom)
                throw boost::program_options::validation_error(boost::program_options::validation_error::invalid_option_value);
        }

        /// Returns string representation of this validation
        static std::string validationString() { return type_string<T>() + u8"<" + (denom == 1 ? output_string(upper) : output_string(upper / (double) denom)); }
};

}}
