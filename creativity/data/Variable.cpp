#include "creativity/data/Variable.hpp"
#include <regex>
#include <sstream>

using namespace Eigen;

namespace creativity { namespace data {

std::string to_string(double d) {
    return std::regex_replace(
            std::regex_replace(std::to_string(d),
                std::regex("(\\.\\d*?)0+$"),
                "$1"),
            std::regex("\\.$"),
            "");
}

VectorXd Variable::values(unsigned int rows, unsigned int offset, unsigned int trim) const {
    if (rows == 0) rows = size();
    VectorXd col(rows);
    populate(col, offset, trim);
    return col;
}

// Returns the unconditionally bracketed name in the base class; the simple classes override below
std::string Variable::nameBracketed(const std::string &left, const std::string &right) const { return left + name() + right; }
std::string SimpleVariable::nameBracketed(const std::string&, const std::string&) const { return name(); }
std::string ConstantVariable::nameBracketed(const std::string&, const std::string&) const { return name(); }

Variable::SizeError::SizeError(const std::string &why) : std::logic_error(why) {}
Variable::SizeError::SizeError()
    : SizeError("populate() size mismatch: target column size + offset + trim != source data size") {}

ConstantVariable::ConstantVariable(double c) : c_{c} {}

void ConstantVariable::populate(Ref<VectorXd> column, unsigned int, unsigned int) const {
    column.setConstant(c_);
}

std::string ConstantVariable::name() const {
    return c_ == 1 ? "const" : c_ == 0 ? "noconst" : to_string(c_);
}

const double& ConstantVariable::value() const { return c_; }

SimpleVariable::SimpleVariable(const std::string &name, const Ref<const VectorXd> values) :
    name_{name}, col_(values)
{}

void SimpleVariable::populate(Ref<VectorXd> column, unsigned int offset, unsigned int trim) const {
    if (column.size() != col_.size() - offset - trim) throw Variable::SizeError();
    column = col_.segment(offset, column.size());
}

std::string SimpleVariable::name() const { return name_; }

unsigned int SimpleVariable::size() const { return col_.size(); }
unsigned int ConstantVariable::size() const { throw SizeError("Cannot call size() on a ConstantVariable"); }


BinaryExpr::BinaryExpr(const std::shared_ptr<const Variable> &left, const std::shared_ptr<const Variable> &right) :
    left_{left}, right_{right} {}

void Multiplication::populate(Eigen::Ref<Eigen::VectorXd> column, unsigned int offset, unsigned int trim) const {
    left_->populate(column, offset, trim);
    column.array() *= right_->values(column.size(), offset, trim).array();
}

std::string Multiplication::name() const {
    // Attempt some bracket collapses where order-of-operations allows
    auto leftptr = left_.get(), rightptr = right_.get();
    std::string left(
            dynamic_cast<const Multiplication*>(leftptr) or
            dynamic_cast<const Division*>(leftptr) or
            dynamic_cast<const Power*>(leftptr) or
            dynamic_cast<const Exponential*>(leftptr) or
            dynamic_cast<const Logarithm*>(leftptr)
            ? left_->name() : left_->nameBracketed());
    std::string op("*");
    std::string right(
            dynamic_cast<const Multiplication*>(rightptr) or
            dynamic_cast<const Division*>(rightptr) or
            dynamic_cast<const Power*>(rightptr) or
            dynamic_cast<const Exponential*>(rightptr) or
            dynamic_cast<const Logarithm*>(rightptr)
            ? right_->name() : right_->nameBracketed());

    if (const ConstantVariable *leftconst = dynamic_cast<const ConstantVariable*>(leftptr)) {
        if (leftconst->value() == -1.0) {
            left.clear(); op = "-";
        }
    }
    if (const ConstantVariable *rightconst = dynamic_cast<const ConstantVariable*>(rightptr)) {
        if (rightconst->value() == -1.0) {
            right = left;
            op = "-";
            left.clear();
        }
    }

    return left + op + right;
}

void Addition::populate(Eigen::Ref<Eigen::VectorXd> column, unsigned int offset, unsigned int trim) const {
    left_->populate(column, offset, trim);
    column.array() += right_->values(column.size(), offset, trim).array();
}

std::string Addition::name() const { return left_->name() + "+" + right_->name(); }

unsigned int BinaryExpr::size() const {
    try { return left_->size(); }
    catch (SizeError&) {
        // If both sides throw, we're a complex operation of just constants (so we want to
        // let this throw):
        return right_->size();
    }
}

UnaryExpr::UnaryExpr(const std::shared_ptr<const Variable> &var) : var_(var) {}
unsigned int UnaryExpr::size() const { return var_->size(); }

void Division::populate(Eigen::Ref<Eigen::VectorXd> column, unsigned int offset, unsigned int trim) const {
    left_->populate(column, offset, trim);
    column.array() /= right_->values(column.size(), offset, trim).array();
}

std::string Division::name() const {
    // This is similar (but not identical) to the Multiplication version: the main
    // difference is in the denominator, which has different precedence rules than the top.

    // Attempt some bracket collapses where order-of-operations allows
    auto leftptr = left_.get(), rightptr = right_.get();
    std::string left(
            dynamic_cast<const Multiplication*>(leftptr) or
            dynamic_cast<const Division*>(leftptr) or
            dynamic_cast<const Power*>(leftptr) or
            dynamic_cast<const Exponential*>(leftptr) or
            dynamic_cast<const Logarithm*>(leftptr)
            ? left_->name() : left_->nameBracketed());
    std::string op("/");
    std::string right(
            dynamic_cast<const Power*>(rightptr) or
            dynamic_cast<const Exponential*>(rightptr) or
            dynamic_cast<const Logarithm*>(rightptr)
            ? right_->name() : right_->nameBracketed());

    if (const ConstantVariable *leftconst = dynamic_cast<const ConstantVariable*>(leftptr)) {
        if (leftconst->value() == -1.0) {
            left.clear(); op = "-";
        }
    }
    if (const ConstantVariable *rightconst = dynamic_cast<const ConstantVariable*>(rightptr)) {
        if (rightconst->value() == -1.0) {
            right = left;
            op = "-";
            left.clear();
        }
    }

    return left + op + right;
}

Power::Power(const std::shared_ptr<const Variable> &var, double power) : UnaryExpr(var), power_{power} {}

void Power::populate(Eigen::Ref<Eigen::VectorXd> column, unsigned int offset, unsigned int trim) const {
    var_->populate(column, offset, trim);
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

std::string Power::name() const { return var_->nameBracketed() + "^" + to_string(power_); }

Exponential::Exponential(double base, const std::shared_ptr<const Variable> &var) : UnaryExpr(var), base_{base} {}

Exponential::Exponential(const std::shared_ptr<const Variable> &var) : Exponential(std::exp(1), var) {}

void Exponential::populate(Eigen::Ref<Eigen::VectorXd> column, unsigned int offset, unsigned int trim) const {
    var_->populate(column, offset, trim);
    if (base_ == std::exp(1))
        column = column.array().exp().matrix();
    else if (base_ == 2)
        column = column.unaryExpr([](double c) { return std::exp2(c); });
    else
        column = column.unaryExpr([this](double c) { return std::pow(base_, c); });
}

std::string Exponential::name() const {
    if (base_ == std::exp(1)) return "exp(" + var_->name() + ")"; // NB: don't need bracketed name here
    else return to_string(base_) + "^" + var_->nameBracketed();
}

std::string Exponential::nameBracketed(const std::string &bracketL, const std::string &bracketR) const {
    if (base_ == std::exp(1)) return "exp(" + var_->name() + ")";
    else return bracketL + to_string(base_) + "^" + var_->nameBracketed() + bracketR;
}

void Logarithm::populate(Eigen::Ref<Eigen::VectorXd> column, unsigned int offset, unsigned int trim) const {
    var_->populate(column, offset, trim);
    column = column.array().log().matrix();
}

std::string Logarithm::name() const {
    return "log(" + var_->name() + ")";
}

std::string Logarithm::nameBracketed(const std::string&, const std::string&) const {
    return name();
}


std::shared_ptr<Multiplication> operator* (const std::shared_ptr<const Variable> &left, const std::shared_ptr<const Variable> &right) {
    return Multiplication::create(left, right);
}

std::shared_ptr<Multiplication> operator* (const std::shared_ptr<const Variable> &v, double c) {
    return Multiplication::create(ConstantVariable::create(c), v);
}
std::shared_ptr<Multiplication> operator* (double c, const std::shared_ptr<const Variable> &v) {
    return Multiplication::create(ConstantVariable::create(c), v);
}
std::shared_ptr<Multiplication> operator- (const std::shared_ptr<const Variable> &v) {
    return -1.0 * v;
}
std::shared_ptr<Addition> operator+ (const std::shared_ptr<const Variable> &left, const std::shared_ptr<const Variable> &right) {
    return Addition::create(left, right);
}
std::shared_ptr<Addition> operator+ (const std::shared_ptr<const Variable> &v, double c) {
    return Addition::create(v, ConstantVariable::create(c));
}
std::shared_ptr<Addition> operator+ (double c, const std::shared_ptr<const Variable> &v) {
    return Addition::create(ConstantVariable::create(c), v);
}
std::shared_ptr<Addition> operator- (const std::shared_ptr<const Variable> &left, const std::shared_ptr<const Variable> &right) {
    return Addition::create(left, -right);
}
std::shared_ptr<Addition> operator- (const std::shared_ptr<const Variable> &v, double c) {
    return Addition::create(v, ConstantVariable::create(-c));
}
std::shared_ptr<Addition> operator- (double c, const std::shared_ptr<const Variable> &v) {
    return c + (-v);
}
std::shared_ptr<Division> operator/ (const std::shared_ptr<const Variable> &numerator, const std::shared_ptr<const Variable> &denominator) {
    return Division::create(numerator, denominator);
}
std::shared_ptr<Division> operator/ (const std::shared_ptr<const Variable> &v, double c) {
    return Division::create(v, ConstantVariable::create(c));
}
std::shared_ptr<Division> operator/ (double c, const std::shared_ptr<const Variable> &v) {
    return Division::create(ConstantVariable::create(c), v);
}
std::shared_ptr<Power> operator^ (const std::shared_ptr<const Variable> &val, double pow) {
    return Power::create(val, pow);
}
std::shared_ptr<Exponential> operator^ (double base, const std::shared_ptr<const Variable> &val) {
    return Exponential::create(base, val);
}

}}

namespace std {
/// std::exp specialization for a Variable.  Returns an Exponential variable wrapper.
shared_ptr<creativity::data::Exponential>
exp(const shared_ptr<const creativity::data::Variable> &var) { return creativity::data::Exponential::create(var); }
/// std::exp2 specialization for a Variable.  Returns an Exponential variable wrapper with base 2.
shared_ptr<creativity::data::Exponential>
exp2(const shared_ptr<const creativity::data::Variable> &var) { return creativity::data::Exponential::create(2, var); }
/// std::log specialization for a Variable.  Returns a Logarithm variable wrapper.
shared_ptr<creativity::data::Logarithm>
log(const shared_ptr<const creativity::data::Variable> &var) { return creativity::data::Logarithm::create(var); }
/// std::sqrt specialization for a Variable.  Returns a Power variable wrapper with power = 0.5.
shared_ptr<creativity::data::Power>
sqrt(const shared_ptr<const creativity::data::Variable> &var) { return creativity::data::Power::create(var, 0.5); }
/// std::pow specialization for a Variable raised to a numeric power.  Returns a Power variable wrapper.
shared_ptr<creativity::data::Power>
pow(const shared_ptr<const creativity::data::Variable> &var, double power) { return creativity::data::Power::create(var, power); }
/// std::pow specialization for a double raised to a Variable.  Returns an Exponential variable wrapper.
shared_ptr<creativity::data::Exponential>
pow(double base, const shared_ptr<const creativity::data::Variable> &var) { return creativity::data::Exponential::create(base, var); }
}
