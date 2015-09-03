#include "creativity/data/Variable.hpp"

using namespace Eigen;

namespace creativity { namespace data {

VectorXd Variable::values(unsigned int rows, unsigned int offset, unsigned int trim) const {
    VectorXd col(rows);
    populate(col, offset, trim);
    return col;
}

Variable::SizeError::SizeError(const std::string &why) : std::logic_error(why) {}
Variable::SizeError::SizeError()
    : SizeError("populate() size mismatch: target column size + offset + trim != source data size") {}

ConstantVariable::ConstantVariable(double constant) : c{constant} {}

void ConstantVariable::populate(Ref<VectorXd> column, unsigned int, unsigned int) const {
    column.setConstant(c);
}

std::string ConstantVariable::name() const {
    return c == 1.0 ? "const" : std::to_string(c);
}

SimpleVariable::SimpleVariable(const std::string &name, const Ref<const VectorXd> &values) :
    name_{name}, col_{values}
{}

void SimpleVariable::populate(Ref<VectorXd> column, unsigned int offset, unsigned int trim) const {
    if (column.size() != col_.size() - offset - trim) throw Variable::SizeError();
    column = col_.segment(offset, column.size());
}

std::string SimpleVariable::name() const { return name_; }


}}
