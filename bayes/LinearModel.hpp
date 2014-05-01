#pragma once
#include <vector>
#include <Eigen/Core>

namespace bayes {

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Dynamic;

/** Class representing a linear regression model with K coefficients (including a constant, if
 * applicable) with Normal-Gamma distribution of coefficients and precision (\f$=
 * \frac{1}{\sigma^2}).
 */
class LinearNG {
    public:
        /** Constructor for a linear model which takes the model parameters.
         *
         * \param beta the vector of beta coefficients.  For a frequentist model, the estimate; for
         * a Bayesian model, the prior mean.
         * \param sigma the standard deviation (or estimate thereof) of the beta estimates.  0 if
         * beta are exact values rather than estimates.
         * \param V the covariance matrix of the `beta` coefficients (or estimate thereof), without .  Note
         * that this estimate is often premultiplied by sigma^2.
         * , *without* the \f$s^2\f$
         * multiple.  In the OLS analog, this is the \f$(X^\top X)^{-1}\f$ matrix of the
         * \f$s^2(X^\top X)^{-1}\f$ covariance estimator.  `V` must be a symmetric,
         * positive-definite matrix, but this is not checked.
         * \param n the size of the prior, typically degrees of freedom.  For example, for an OLS
         * prior with 50 observations, this should be 50-K.
         */
        LinearModel(VectorXd beta, double sigma, MatrixXd V, double n)
            beta_{beta}, 
            : prior_beta_{std::move(beta)}, prior_s_sq_{std::move(s_sq)}, prior_v_{std::move(V)}, prior_n_{std::move(n)}
        {}

        /** Constructs a LinearModel which uses the posterior of another LinearModel object.
         */
        //static LinearModel fromPosterior(const LinearModel &lm);

        /// Accesses the prior beta values
        const Matrix<double, K, 1>& priorBeta() const { return prior_beta_; }
        /// Accesses the prior s^2 value
        const double& priorSSq() const { return prior_s_sq_; }
        /// Accesses the prior V value
        const Matrix<double, K, K>& priorV() const { return prior_v_; }
        /// Access the prior size
        const double& priorN() const { return prior_n_; }

        /** Access the matrix of beta draws, not including burn draws.  Each column is a set of K
         * beta values.
         */
        const Matrix<double, K, Dynamic>& betaDraws() { return draw_beta_; }
        /** Access the row vector of s^2 values, not including burn draws.
         */
        const VectorXd& sSqDraws() { return draw_s_sq_; }


    private:
        /// Priors:
        const Matrix<double, K, 1> prior_beta_;
        const double prior_s_sq_;
        const Matrix<double, K, K> prior_v_;
        const double prior_n_;

        /// Draws.  These begin life as null matrices.
        Matrix<double, K, Dynamic> draw_beta_{};
        VectorXd draw_s_sq_{};

        template <size_t K2, size_t N, class R>
        friend class Gibbs;

};

}
