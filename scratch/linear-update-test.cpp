#include "creativity/belief/Linear.hpp"
#include <eris/Random.hpp>
#include <Eigen/Core>
#include <iostream>

#define print_model(M) \
    std::cout << #M ":\n" << \
        "beta_: " << M.beta().transpose() << "\n" << \
        "n_: " << M.n() << "\n" << \
        "s2_: " << M.s2() << "\n" << \
        "V_:\n" << M.V() << "\n" << \
        "\n"

using namespace creativity::belief;
using namespace Eigen;
using namespace eris;

int main() {

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

    Linear foo_10_oneshot = foo.update(y.head(10), X.topRows(10));

    Linear foo_10_twoshot = foo.update(y.head(5), X.topRows(5));
    foo_10_twoshot = foo_10_twoshot.update(y.middleRows(5, 5), X.middleRows(5, 5));

    Linear foo_10_tenshot = foo;
    for (int i = 0; i < 10; i++) {
        foo_10_tenshot = foo_10_tenshot.update(y.row(i), X.row(i));
    }

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

    print_model(foo_10_oneshot);
    print_model(foo_10_twoshot);
    print_model(foo_10_tenshot);

    print_model(foo_100_oneshot);
    print_model(foo_100_twoshot);
    print_model(foo_100_fiveshot);
    print_model(foo_100_tenshot);
    print_model(foo_100_hundredshot);
}
