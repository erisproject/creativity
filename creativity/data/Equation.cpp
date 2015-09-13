#include "creativity/data/Equation.hpp"

namespace creativity { namespace data {

Equation::Proxy::Proxy(Equation &eq) : eq_(eq) {}

Equation::Equation(const Variable &y) : dep_var_(y) {}

Equation::Proxy& Equation::Proxy::operator+(const Variable &var) { eq_.addVar(var); return *this; }
Equation::Proxy Equation::operator%(const Variable &var) { addVar(var); return Proxy(*this); }
Equation::Proxy Equation::operator%(double c) { return *this % ConstantVariable(c); }

Equation Equation::operator+(const Variable &var) const & { Equation copy(*this); copy.addVar(var); return copy; }
Equation Equation::operator+(const Variable &var) && { addVar(var); return std::move(*this); }
Equation Equation::operator+(double c) const & { return *this + ConstantVariable(c); }
Equation Equation::operator+(double c) && { return std::move(*this) + ConstantVariable(c); }

void Equation::addVar(const Variable &var) {
    indep_vars_.emplace_back(var);
}

template <>
void Equation::addVar<ConstantVariable>(ConstantVariable &&v) {
    const_ = std::move(v);
}

const Variable& Equation::depVar() const { return dep_var_; }

std::vector<std::string> Equation::names() const {
    std::vector<std::string> names;
    for (const Variable &v : *this) {
        names.push_back(v.name());
    }
    return names;
}

unsigned int Equation::numVars() const {
    return hasConstant()
        ? indep_vars_.size()
        : indep_vars_.size() - 1;
}

bool Equation::hasConstant() const {
    const ConstantVariable &constant = static_cast<const ConstantVariable&>((const Variable&) indep_vars_.front());
    return constant.value() != 0;
}

const std::list<std::reference_wrapper<const Variable>>::const_iterator Equation::begin() const {
    auto it = indep_vars_.begin();
    if (not hasConstant()) it++;
    return it;
}

/** Const access to the past-the-end iterator of independent variables.
 */
const std::list<std::reference_wrapper<const Variable>>::const_iterator Equation::end() const {
    return indep_vars_.end();
}

std::ostream& operator<<(std::ostream &os, const Equation &eq) {
    os << eq.dep_var_.nameBracketed("[", "]");
    bool first = true;
    for (const Variable &var : eq.indep_vars_) {
        if (first) { os << " ~ "; first = false; }
        else { os << " + "; }
        os << var.nameBracketed("[", "]");
    }
    return os;
}

}}
