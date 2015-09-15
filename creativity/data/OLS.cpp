#include "creativity/data/OLS.hpp"
#include "creativity/data/Variable.hpp"
#include "creativity/data/Regression.hpp"
#include "creativity/data/tabulate.hpp"
#include <Eigen/SVD>
#include <algorithm>
#include <limits>
#include <string>
#include <ostream>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/complement.hpp>

using namespace Eigen;
namespace creativity { namespace data {

OLS::OLS(const Equation &model) : model_(model) {}
OLS::OLS(Equation &&model) : model_(std::move(model)) {}

void OLS::gather() {
    if (gathered_) return;

    y_ = model_.depVar().values();
    X_ = MatrixXd(y_.size(), model_.numVars());

    unsigned int k = 0;
    for (const Variable &var : model_) {
        var.populate(X_.col(k++));
    }

    gathered_ = true;
}

const MatrixXd& OLS::X() { gather(); return X_; }
const VectorXd& OLS::y() { gather(); return y_; }

void OLS::solve() {
    if (solved_) return;

    gather();
    const unsigned int n = X_.rows(), k = X_.cols();

    if (n < k) throw RankError("Cannot compute OLS estimates with n < k");

    JacobiSVD<MatrixXd> svd(X_, ComputeThinU | ComputeThinV);
    if (svd.rank() < X_.cols()) throw RankError("Cannot compute OLS estimates: independent data is not full column rank (" + std::to_string(svd.rank()) + " < " + std::to_string(X_.cols()) + ")");

    beta_ = svd.solve(y_);

    residuals_ = y_ - X_ * beta_;

    ssr_ = residuals_.squaredNorm();
    s2_ = X_.cols() == X_.rows() ? std::numeric_limits<double>::quiet_NaN() : ssr_ / (n - k);
    if (model_.hasConstant()) {
        // Centered R^2:
        double sum_y = y_.sum();
        R2_ = 1 - ssr_ / (y_.squaredNorm() - sum_y*sum_y / n);
    }
    else {
        // Uncentered R^2
        R2_ = 1 - ssr_ / y_.squaredNorm();
    }

    DiagonalMatrix<double, Dynamic> S2i(k);
    // We're after: V S^{-2} V', where S is the (diagonal) matrix of SVD singular values and V is
    // the SVD "V" matrix.  Computing the ^{-2} on S is trivial, since it's diagonal: just compute
    // ^{-2} of each element.
    S2i.diagonal() = svd.singularValues().array().square().inverse();
    var_beta_ = s2_ * svd.matrixV() * S2i * svd.matrixV().transpose();

    se_ = var_beta_.diagonal().cwiseSqrt();
    t_ratios_ = beta_.array() / se_.array();
    boost::math::students_t tdist(n - k);
    p_values_ = t_ratios_.unaryExpr([&tdist](double t) { return 2.0*cdf(tdist, t < 0 ? t : -t); });

    solved_ = true;
}

const VectorXd& OLS::beta() const { requireSolved(); return beta_; }
const MatrixXd& OLS::covariance() const { requireSolved(); return var_beta_; }
const VectorXd& OLS::residuals() const { requireSolved(); return residuals_; }
const double& OLS::s2() const { requireSolved(); return s2_; }
const double& OLS::ssr() const { requireSolved(); return ssr_; }
const double& OLS::Rsq() const { requireSolved(); return R2_; }
const VectorXd& OLS::se() const { requireSolved(); return se_; }
const VectorXd& OLS::tRatios() const { requireSolved(); return t_ratios_; }
const VectorXd& OLS::pValues() const { requireSolved(); return p_values_; }

OLS::ftest OLS::fTest() const {
    requireSolved();
    ftest result;
    result.df_numerator = k();
    result.df_denominator = df();
    double rssr = y_.squaredNorm();
    if (model_.hasConstant()) {
        result.df_numerator--;
        double sum_y = y_.sum();
        rssr -= sum_y * sum_y / n();
    }
    result.f = (rssr - ssr_) * result.df_denominator / (ssr_ * result.df_numerator);
    boost::math::fisher_f fdist(result.df_numerator, result.df_denominator);
    result.p = cdf(complement(fdist, result.f));
    return result;
}


std::ostream& operator<<(std::ostream &os, const OLS &ols) {
    os << "OLS for model: " << ols.model();

    if (ols.solved_) {
        Eigen::MatrixXd results(ols.beta_.size(), 4);
        results.col(0) = ols.beta();
        results.col(1) = ols.se();
        results.col(2) = ols.tRatios();
        results.col(3) = ols.pValues();

        std::vector<std::string> rownames;
        for (const Variable &var : ols.model()) {
            rownames.emplace_back(var.name());
        }

        std::vector<std::string> stars;
        stars.push_back("");
        for (int i = 0; i < ols.pValues().size(); i++) {
            const double &p = ols.pValues()[i];
            stars.push_back(
                    p < .001 ? "***" :
                    p < .01 ? "**" :
                    p < .05 ? "*" :
                    p < .1 ? "." :
                    "");
        }

        os << "; n=" << ols.n() << " observations\n" << tabulate(results, TableFormat::Text, rownames, {"Coefficient", "std.err.", "t-stat", "p-value"}, stars);
        os << "---\nSignificance: 0 `***' .001 `**' .01 `*' .05 `.' .1 ` ' 1\n\n";
        os << "Residual s.e.: " << std::sqrt(ols.s2()) << " on " << ols.df() << " d.f.\n";
        auto ftest = ols.fTest();
        os << "R^2: " << ols.Rsq() << "\n";
        os << "F-statistic: " << ftest.f << " on " << ftest.df_numerator << " and " << ftest.df_denominator << " d.f., p-value = " << ftest.p << "\n";

    }
    else {
        os << " (not solved)";
    }
    return os;
}

}}
