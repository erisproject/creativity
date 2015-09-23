#include "creativity/data/Equation.hpp"

#include <iostream>

namespace creativity { namespace data {

Equation::Proxy::Proxy(Equation &eq) : eq_(eq) {}

Equation::Equation(const std::shared_ptr<const Variable> &y) : dep_var_(y) {}

Equation::Proxy& Equation::Proxy::operator+(const std::shared_ptr<const Variable> &var) { eq_.addVar(var); return *this; }
Equation::Proxy Equation::operator%(const std::shared_ptr<const Variable> &var) { addVar(var); return Proxy(*this); }
Equation::Proxy Equation::operator%(double c) { return *this % ConstantVariable::create(c); }

Equation Equation::operator+(const std::shared_ptr<const Variable> &var) const & { Equation copy(*this); copy.addVar(var); return copy; }
Equation Equation::operator+(const std::shared_ptr<const Variable> &var) && { addVar(var); return std::move(*this); }
Equation Equation::operator+(double c) const & { return *this + ConstantVariable::create(c); }
Equation Equation::operator+(double c) && { return std::move(*this) + ConstantVariable::create(c); }

void Equation::addVar(const std::shared_ptr<const Variable> &var) {
    if (dynamic_cast<const ConstantVariable*>(var.get()))
        indep_vars_.front() = var;
    else
        indep_vars_.push_back(var);
}

std::shared_ptr<const Variable> Equation::depVar() const { return dep_var_; }

std::vector<std::string> Equation::names() const {
    std::vector<std::string> names;
    for (const auto &v : *this) {
        names.push_back(v->name());
    }
    return names;
}

unsigned int Equation::numVars() const {
    return hasConstant()
        ? indep_vars_.size()
        : indep_vars_.size() - 1;
}

bool Equation::hasConstant() const {
    return dynamic_cast<const ConstantVariable&>(*(indep_vars_.front())).value() != 0;
}

std::list<std::shared_ptr<const Variable>>::const_iterator Equation::begin() const {
    auto it = indep_vars_.begin();
    if (not hasConstant()) it++;
    return it;
}

/** Const access to the past-the-end iterator of independent variables.
 */
std::list<std::shared_ptr<const Variable>>::const_iterator Equation::end() const {
    return indep_vars_.end();
}

std::ostream& operator<<(std::ostream &os, const Equation &eq) {
    os << eq.dep_var_->nameBracketed("[", "]");
    bool first = true;
    for (const auto &var : eq.indep_vars_) {
        if (first) { os << " ~ "; first = false; }
        else { os << " + "; }
        os << std::flush;
        os << var->nameBracketed("[", "]");
        os << std::flush;
    }
    return os;
}

}}
