#include "creativity/data/Variable.hpp"
#include "creativity/data/Equation.hpp"
#include "creativity/data/OLS.hpp"

#include <Eigen/Core>
#include <iostream>
#define ERIS_DEBUG
#include <eris/debug.hpp>

using namespace creativity::data;
using namespace Eigen;

int main(int, char **) {

    Eigen::VectorXd ydata(5);
    ydata << 1, 4, 4, 3.2, -1.7;
    Eigen::MatrixXd Xdata(5, 3);
    Xdata << 1, 2, 3,
        4, 4, 4,
        -5, 8, -5,
        100, 10, 1.7,
        0, 1e-10, -8.4;

    std::cout << Xdata << "\n";

    SimpleVariable y("y", ydata);
    SimpleVariable c1("X_1", Xdata.col(0));
    SimpleVariable c2("X_2", Xdata.col(1));
    SimpleVariable c3("X_3", Xdata.col(2));
    ConstantVariable constant(1);

    auto c1sq = c1 ^ 2;
    auto c2exp = std::exp(c2);
    auto const2 = constant * 2;

    ERIS_DBGVAR(c1sq.values(4, 0, 1));
    ERIS_DBGVAR(c2exp.values(4, 1));
    ERIS_DBGVAR(const2.values(7));

    ERIS_DBGVAR((c1 * c3).values(5));
    ERIS_DBGVAR((c1 + c3).values(5));
    ERIS_DBGVAR((c1 - c3).values(5));
    ERIS_DBGVAR((c1 - 3).values(5));
    ERIS_DBGVAR(-c2.values(5));

    Equation model1(SimpleVariable("y", ydata));
    model1 % 1 + c1 + c2;

    Equation model2(y);
    model2 % 0; model2 % c1; model2 % c2 + c3;

    Equation model3(y);

    Equation model4(y);
    model4 % 2 + 2*c1 + (2*c2 + 2*c3);

    model1 % std::exp(c3) + std::exp2(c2) + std::log(c2) + std::pow(3, std::log(std::exp(c2)));

    std::cout << "Model1: " << model1 << "\nModel2: " << model2 << "\nModel3: " << model3 << "\nModel4: " << model4 << "\n";

    std::cout << "Running model4:\n";
    OLS ols4(model4);
    ols4.solve();
    std::cout << ols4 << "\n";
}
