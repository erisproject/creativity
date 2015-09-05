#include "creativity/data/Variable.hpp"
#include <regex>

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

SimpleVariable::SimpleVariable(const std::string &name, const Ref<const VectorXd> &values) :
    name_{name}, col_{values}
{}

void SimpleVariable::populate(Ref<VectorXd> column, unsigned int offset, unsigned int trim) const {
    if (column.size() != col_.size() - offset - trim) throw Variable::SizeError();
    column = col_.segment(offset, column.size());
}

std::string SimpleVariable::name() const { return name_; }

unsigned int SimpleVariable::size() const { return col_.size(); }
unsigned int ConstantVariable::size() const { throw SizeError("Cannot call size() on a ConstantVariable"); }

}}
