#include "bayes/PosteriorSimulator.hpp"

namespace bayes {


PosteriorSimulator::PosteriorSimulator(const MatrixXd &Xa, const MatrixXd &ya)
    : X{Xa}, y{ya}
{
    if (X.rows() != y.rows())
        throw std::runtime_error("Error in PosteriorSimulator(X,y): X.rows != y.rows");
}

PosteriorSimulator::PosteriorSimulator(const MatrixXd &Xa) : X{Xa}, y{Eigen::Matrix<double,0,0>()}
{}


void PosteriorSimulator::burn(unsigned int b) {
    for (unsigned int i = 0; i < b; i++)
        discard();
}

void PosteriorSimulator::discard() {
    draw();
}

MatrixXd PosteriorSimulator::drawMany(unsigned int s) {
    MatrixXd draws;
    for (unsigned int i = 0; i < s; i++) {
        draws.row(i) = draw();
    }
    return draws;
}

}
