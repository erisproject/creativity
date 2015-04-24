#include "creativity/belief/Linear.hpp"
#include <eris/Random.hpp>
#include <Eigen/Core>
#include <Eigen/QR>
#include <iostream>
#include <iomanip>

#define print_model(M) \
    std::cout << #M ":\n" << \
        "beta_: " << M.beta().transpose() << "\n" << \
        "n_: " << M.n() << "\n" << \
        "s2_: " << M.s2() << "\n" << \
        "V_inv_:\n" << M.Vinv() << "\n" << \
        "\n"

using namespace creativity::belief;
using namespace Eigen;
using namespace eris;

int main() {

    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10) << std::fixed;
    Linear foo(3);

    VectorXd beta(3);
    beta << -1, 4, 0.5;
    MatrixXd X(100, 3);
    VectorXd y(100);
    VectorXd u(100);

    auto &rng = Random::rng();
    std::normal_distribution<double> stdnormal;

    for (int t = 0; t < 100; t++) {
        for (int k = 0; k < 3; k++) {
            X(t,k) = stdnormal(rng);
        }
        u[t] = 2.5*stdnormal(rng);
    }

    y = X * beta + u;

    VectorXd betahat = (X.transpose() * X).fullPivHouseholderQr().solve(X.transpose() * y);
    std::cout << "OLS:\n" << "beta^: " << betahat.transpose() << "\n" <<
        "n: " << X.rows() << "\n" <<
        "sigmahat^2: " << (y - X * betahat).squaredNorm() / X.rows() << "\n" <<
        "X'X:\n" << (X.transpose() * X) << "\n" <<
        "\n";

    Linear foo_100_oneshot = foo.update(y, X);

    Linear foo_100_twoshot = foo.update(y.topRows(50), X.topRows(50));
    foo_100_twoshot = foo_100_twoshot.update(y.bottomRows(50), X.bottomRows(50));

    Linear foo_100_fiveshot = foo;
    for (int i = 0; i < 100; i += 20) {
        foo_100_fiveshot = foo_100_fiveshot.update(y.middleRows(i, 20), X.middleRows(i, 20));
    }

    Linear foo_100_tenshot = foo;
    for (int i = 0; i < 10; i++) {
        foo_100_tenshot = foo_100_tenshot.update(y.middleRows(10*i, 10), X.middleRows(10*i, 10));
    }

    Linear foo_100_hundredshot = foo;
    for (int i = 0; i < 100; i++) {
        foo_100_hundredshot = foo_100_hundredshot.update(y.row(i), X.row(i));
    }

    print_model(foo_100_oneshot);
    print_model(foo_100_twoshot);
    print_model(foo_100_fiveshot);
    print_model(foo_100_tenshot);
    print_model(foo_100_hundredshot);

    Linear foo_100_weakened_fiftyshot = foo.update(y.head(50), X.topRows(50));
    foo_100_weakened_fiftyshot = foo_100_weakened_fiftyshot.weaken(2).update(y.tail(50), X.bottomRows(50));

    MatrixXd Xw = X;
    Xw.topRows(50) *= 0.5;
    VectorXd yw = y;
    yw.head(50) *= 0.5;
    Linear foo_100_weakened_direct = foo.update(yw, Xw);

    print_model(foo_100_weakened_direct);
    print_model(foo_100_weakened_fiftyshot);
}
