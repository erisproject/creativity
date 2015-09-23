#include "creativity/data/Variable.hpp"
#include "creativity/data/Equation.hpp"
#include "creativity/data/OLS.hpp"
#include "creativity/data/SUR.hpp"

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

    auto y = SimpleVariable::create("y", ydata);
    auto c1 = SimpleVariable::create("X_1", Xdata.col(0).transpose());
    auto c2 = SimpleVariable::create("X_2", Xdata.col(1).transpose());
    auto c3 = SimpleVariable::create("X_3", Xdata.col(2).transpose());
    auto constant = ConstantVariable::create(1);

    auto c1sq = c1 ^ 2;
    auto c2exp = std::exp(c2);
    auto const2 = constant * 2;

    ERIS_DBGVAR(c1sq->values(4, 0, 1).transpose());
    ERIS_DBGVAR(c2exp->values(4, 1).transpose());
    ERIS_DBGVAR(const2->values(7).transpose());

    ERIS_DBGVAR((c1 * c3)->values(5).transpose());
    ERIS_DBGVAR((c1 + c3)->values(5).transpose());
    ERIS_DBGVAR((c1 - c3)->values(5).transpose());
    ERIS_DBGVAR((c1 - 3)->values(5).transpose());
    ERIS_DBGVAR(((-c2)->values(5)).transpose());

    Equation model1(SimpleVariable::create("y", ydata));
    model1 % 1 + c1 + c2;

    Equation model2(y);
    model2 % 0; model2 % c1; model2 % c2 + c3;

    Equation model3(y);

    Equation model4(y);
    model4 % 2 + 2*c1 + (2*c2 + 2*c3);

    model1 % std::exp(c3) + std::exp2(c2) + std::log(c2) + std::pow(3, std::log(std::exp(c2)));

    std::cout << "Model1: " << model1 << "\nModel2: " << model2 << "\nModel3: " << model3 << "\nModel4: " << model4 << "\n";

    OLS ols4(model4);
    std::cout << "ols4, pre-solve:\n" << ols4 << "\n";
    ols4.solve();
    std::cout << "ols4, solved:\n" << ols4 << "\n";

    VectorXd testbeta(4), testu(10), testy(10);
    MatrixXd testX(10, 4);
    testbeta << -3, 1, -10, 0.5;
    testu << 1.48683859, -1.17687648, -0.68941946, -0.30323162, -0.03485257,  3.57557769, -0.07831670, -1.22086850, -1.65255339,  2.11408666;
    testX <<
        1, -0.6820301, -0.713643322, -0.9163320,
        1,  0.2300387,  0.870683685, -1.4051567,
        1,  0.3266254, -0.191723819,  0.3227581,
        1,  0.3741961, -0.304324581, -0.6928428,
        1, -2.4596900,  0.376452159,  2.8230870,
        1, -0.3383137, -0.007462838,  0.1665202,
        1, -1.3131063, -0.192922778, -0.3996082,
        1,  0.4725750,  0.681960046,  0.4801447,
        1,  0.3315413, -0.411196959, -0.9939557,
        1, -0.4342752,  0.778097556,  0.6067932;
    testy = testX * testbeta + testu;
    OLS testols(Equation(SimpleVariable::create("y", testy))
            + SimpleVariable::create("X1", testX.col(1))
            + SimpleVariable::create("X2", testX.col(2))
            + SimpleVariable::create("X3", testX.col(3)));
    testols.solve();
    std::cout << "Another OLS test:\n" << testols << "\n";

    SUR sur4(model4);
    std::cout << "sur4, pre-solve:\n" << sur4 << "\n";
    sur4.solve();
    std::cout << "sur4, solved:\n" << sur4 << "\n";

    SUR sur9(model4, Equation(c1) + 0 + c2);
    sur9.solve();
    std::cout << sur9 << "\n";
}
