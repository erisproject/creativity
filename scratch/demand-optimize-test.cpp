#include "creativity/belief/Demand.hpp"
#include <iostream>
#include <iomanip>

using Eigen::VectorXd;
using Eigen::MatrixXd;

int main() {

    VectorXd beta{7};
    beta << 15, -0.05, 0.17, -3, -0.0246, 0.112, -0.123;

    MatrixXd V = MatrixXd::Zero(7, 7); // A perfectly informed prior.  Nonsense, of course, but I'm not testing Bayesian updating here.

    creativity::belief::Demand demand{3, beta, 0, V, 1e10};

    double q = 13;
    unsigned long S = 40, other = 3, market = 5;
    
    std::cout << std::setprecision(16);
    std::cout << "optimal for c=0: " << demand.argmaxP(q, S, other, market) << "\n";
    std::cout << "optimal for c=1e-100: " << demand.argmaxP(q, S, other, market, 1e-100) << "\n";
}

