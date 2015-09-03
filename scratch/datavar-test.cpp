#include "creativity/data/Variable.hpp"
#include <Eigen/Core>
#include <iostream>
#define ERIS_DEBUG
#include <eris/debug.hpp>

using namespace creativity::data;
using namespace Eigen;

int main(int, char **) {

    Eigen::MatrixXd foo(5, 3);
    foo << 1, 2, 3,
        4, 4, 4,
        -5, 8, -5,
        10, 1e100, std::numeric_limits<double>::quiet_NaN(),
        0, 0, std::numeric_limits<double>::infinity();

    std::cout << foo << "\n";

    SimpleVariable c1("X_1", foo.col(0));
    SimpleVariable c2("X_2", foo.col(1));
    SimpleVariable c3("X_3", foo.col(2));
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

}
