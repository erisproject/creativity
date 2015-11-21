#include "creativity/data/SUR.hpp"
#include "creativity/data/Regression.hpp"
#include "creativity/data/Variable.hpp"
#include "creativity/data/tabulate.hpp"
#include <Eigen/SVD>
#include <Eigen/QR>
#include <stdexcept>
#include <string>
#include <ostream>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/complement.hpp>

using namespace Eigen;
namespace creativity { namespace data {

const std::vector<Equation>& SUR::equations() const { return eqs_; }

void SUR::clear() {
    gathered_ = false;
    X_.clear();
    y_.resize(0);
    solved_ = false;
    beta_.resize(0);
    var_beta_.resize(0, 0);
    residuals_.resize(0);
    ssr_.clear();
    s2_.clear();
    R2_.clear();
    t_ratios_.resize(0);
    p_values_.resize(0);
    // FIXME: other vars?
}

void SUR::gather() {
    if (gathered_) return;

    if (eqs_.empty()) throw std::logic_error("SUR::gather(): Invalid SUR: no equations specified");

    // Since this is an SUR model, each model should have the same number of observations, so check
    // that.  The number of regressors, on the other hand, is the sum of the number of regressors in
    // each equation (which can differ between equations).
    unsigned int per_eq_rows = eqs_.front().depVar()->size();
    unsigned int total_cols = 0;
    for (const auto& eq : eqs_) {
        if (eq.depVar()->size() != per_eq_rows) throw Variable::SizeError("Invalid SUR: found equations with differing numbers of rows");
        total_cols += eq.numVars();
    }
    unsigned int total_rows = per_eq_rows * eqs_.size();

    y_.resize(total_rows);
    X_.reserve(eqs_.size());
    // Second pass: actually populate
    unsigned int r = 0, total_c = 0;
    for (const auto& eq : eqs_) {
        eq.depVar()->populate(y_.segment(r, per_eq_rows));
        X_.emplace_back(per_eq_rows, eq.numVars());
        int c = 0;
        for (const auto &var : eq) {
            var->populate(X_.back().col(c++));
        }
        total_c += c;
        r += per_eq_rows;
    }
    if (r != total_rows or total_c != total_cols) {
        throw std::runtime_error("Internal error: unexpected number of rows/columns");
    }

    gathered_ = true;
}

const unsigned& SUR::k(unsigned i) const { return k_[i]; }
int SUR::df(unsigned i) const { return n() - k(i); }

const std::vector<MatrixXd>& SUR::X() const { requireGathered(); return X_; }
const VectorXd& SUR::y() const { requireGathered(); return y_; }
VectorBlock<const VectorXd> SUR::y(unsigned i) const {
    requireGathered();
    return VectorBlock<const VectorXd>(y_, i*n(), n());
}

void SUR::solve() {
    if (solved_) return;

    gather();

    // Step 1: run OLS
    const unsigned per_eq_rows = n(), num_eqs = eqs_.size(), K = offset_.back() + k_.back();

    // FIXME: might want to iterate on the residuals!

    // First step: estimate each equation by OLS to get OLS residuals; stack these
    // num_eq residual vectors into a matrix
    MatrixXd residuals(per_eq_rows, num_eqs);
    for (unsigned g = 0; g < num_eqs; g++) {
        JacobiSVD<MatrixXd> svd(X_[g], ComputeThinU | ComputeThinV);
#if EIGEN_VERSION_AT_LEAST(3,2,2)
// svd.rank() was added in 3.2.2; for previous versions just skip the rank safety test
        if (svd.rank() < k_[g]) throw RankError("Cannot compute SUR estimates: equation " + std::to_string(g) + " independent data is not full column rank (" + std::to_string(svd.rank()) + " < " + std::to_string(k_[g]) + ")");
#endif

        auto y_g = y(g);
        residuals.col(g) = y_g - X_[g] * svd.solve(y_g);
    }

    // Now get the Sigma estimate from the reshaped residuals, and its inverse
    MatrixXd Sigmahatinv = (residuals.transpose() * residuals / per_eq_rows).fullPivHouseholderQr().inverse();

    MatrixXd sXtX(K, K);
    VectorXd sXty = VectorXd::Zero(K);
    // c.f. ETM equation (12.14)
    for (unsigned i = 0; i < num_eqs; i++) {
        for (unsigned j = 0; j < num_eqs; j++) {
            // Set up per-inverted, first KxK matrix in (12.14)
            sXtX.block(offset_[i], offset_[j], k_[i], k_[j]) = Sigmahatinv(i,j) * (X_[i].transpose() * X_[j]);

            if (i != j) sXtX.block(offset_[j], offset_[i], k_[j], k_[i]) =
                sXtX.block(offset_[i], offset_[j], k_[i], k_[j]).transpose();
        }
        for (unsigned j = 0; j < num_eqs; j++) {
            // Set up right-hand vector expression in (12.14)
            sXty.segment(offset_[i], k_[i]) += Sigmahatinv(i,j) * X_[i].transpose() * y_.segment(j*per_eq_rows, per_eq_rows);
        }
    }

    var_beta_ = sXtX.fullPivHouseholderQr().inverse();
    beta_ = var_beta_ * sXty;
    se_= var_beta_.diagonal().cwiseSqrt();

    residuals_.resize(per_eq_rows * num_eqs);
    for (unsigned i = 0; i < num_eqs; i++) {
        residuals_.segment(i*per_eq_rows, per_eq_rows) = y(i) - X_[i] * beta_.segment(offset_[i], k_[i]);
    }
    ssr_.resize(num_eqs, std::numeric_limits<double>::quiet_NaN());
    s2_.resize(num_eqs, std::numeric_limits<double>::quiet_NaN());
    R2_.resize(num_eqs, std::numeric_limits<double>::quiet_NaN());
    t_ratios_ = beta_.array() / se_.array();
    p_values_.resize(beta_.size());
    for (unsigned j = 0; j < num_eqs; j++) {
        unsigned start = j*per_eq_rows;
        ssr_[j] = residuals_.segment(start, per_eq_rows).squaredNorm();
        s2_[j] = ssr_[j] / per_eq_rows;
        auto &eq = eqs_[j];
        if (eq.hasConstant()) {
            // Centered R^2:
            double sum_y = y_.segment(start, per_eq_rows).sum();
            R2_[j] = 1 - ssr_[j] / (y_.segment(start, per_eq_rows).squaredNorm() - sum_y*sum_y / per_eq_rows);
        }
        else {
            // Uncentered R^2
            R2_[j] = 1 - ssr_[j] / y_.segment(start, per_eq_rows).squaredNorm();
        }

        boost::math::students_t tdist(per_eq_rows - k_[j]);
        p_values_.segment(offset_[j], k_[j]) = t_ratios_.segment(offset_[j], k_[j]).unaryExpr(
                [&tdist](double t) { return std::isnan(t) ? t : 2.0*cdf(tdist, t < 0 ? t : -t); });
    }

    solved_ = true;
}

VectorBlock<const VectorXd> SUR::beta(unsigned i) const {
    requireSolved();
    return VectorBlock<const VectorXd>(beta_, offset_[i], k_[i]);
}
Block<const MatrixXd> SUR::covariance(unsigned i) const {
    requireSolved();
    return Block<const MatrixXd>(var_beta_, offset_[i], offset_[i], k_[i], k_[i]);
}
VectorBlock<const VectorXd> SUR::residuals(unsigned i) const {
    requireSolved();
    return VectorBlock<const VectorXd>(residuals_, i*n(), n());
}

const double& SUR::s2(unsigned i) const { requireSolved(); return s2_[i]; }
const double& SUR::ssr(unsigned i) const { requireSolved(); return ssr_[i]; }
const double& SUR::Rsq(unsigned i) const { requireSolved(); return R2_[i]; }
VectorBlock<const VectorXd> SUR::se(unsigned i) const {
    requireSolved();
    return VectorBlock<const VectorXd>(se_, offset_[i], k_[i]);
}
VectorBlock<const VectorXd> SUR::tRatios(unsigned i) const {
    requireSolved();
    return VectorBlock<const VectorXd>(t_ratios_, offset_[i], k_[i]);
}
VectorBlock<const VectorXd> SUR::pValues(unsigned i) const {
    requireSolved();
    return VectorBlock<const VectorXd>(p_values_, offset_[i], k_[i]);
}

std::ostream& operator<<(std::ostream &os, const SUR &sur) {
    os << "SUR system; " << sur.equations().size() << " equations, n = " << sur.n() << " observations per eq.";
    if (not sur.solved_) {
        os << " (not solved):\n";
        for (auto &eq : sur.equations()) {
            os << eq << "\n";
        }
    }
    else {
        os << ":\n\n";
        for (unsigned j = 0; j < sur.eqs_.size(); j++) {
            auto &eq = sur.equations()[j];
            std::ostringstream title;
            title << "Equation " << j+1 << ": " << eq << ":";
            Eigen::MatrixXd results(sur.k(j), 4);
            results.col(0) = sur.beta(j);
            results.col(1) = sur.se(j);
            results.col(2) = sur.tRatios(j);
            results.col(3) = sur.pValues(j);

            std::vector<std::string> rownames;
            for (const auto &var : sur.eqs_[j]) {
                rownames.emplace_back(var->name());
            }

            std::vector<std::string> stars;
            stars.push_back("");
            auto pvalues = sur.pValues(j);
            for (int i = 0; i < pvalues.size(); i++) {
                const double &p = pvalues[i];
                stars.push_back(
                        p < .001 ? " ***" :
                        p < .01 ? " **" :
                        p < .05 ? " *" :
                        p < .1 ? " ." :
                        " ");
            }

            os << tabulate(results, tabulation_options(title.str(), TableFormat::Text, os.precision(), "\t"), rownames,
                    {"Coefficient", "std.err.", "t-stat", "p-value"}, stars);

            double mean_y = sur.y(j).mean();
            double sd_y = std::sqrt((sur.y(j).squaredNorm() - sur.n() * mean_y * mean_y) / (sur.n()-1));
            double ssr = sur.ssr(j);
            double sereg = std::sqrt(sur.s2(j));
            double R2 = sur.Rsq(j);

            size_t w1 = 0, w2 = 0;
            for (auto &d : {mean_y, ssr, R2}) {
                std::ostringstream out;
                out.precision(os.precision());
                out << d;
                w1 = std::max(w1, out.str().length());
            }
            for (auto &d : {sd_y, sereg}) {
                std::ostringstream out;
                out.precision(os.precision());
                out << d;
                w2 = std::max(w2, out.str().length());
            }
            os << "\n";
            os << "\tMean dependent var  " << std::setw(w1) << mean_y << "  " << "S.D. dependent var  " << std::setw(w2) << sd_y  << "\n";
            os << "\tSum squared resid   " << std::setw(w1) << ssr    << "  " << "S.E. of regression  " << std::setw(w2) << sereg << "\n";
            os << "\tR-squared           " << std::setw(w1) << R2 << "\n";
            os << "\n";
        }
    }

    return os;
}

}}
