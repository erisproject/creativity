#pragma once
#include <bayes/LinearModel.hpp>
#include <Eigen/LU>
#include <Eigen/Cholesky>
#include <random>

namespace bayes {

using Eigen::Matrix;
using Eigen::VectorXd;

class Gibbs : public PosteriorSimulator {
    public:

        /** Constructs a new Gibbs sampler for the given LinearModel.
         *
         * \param lm a freshly-constructed (or moved) LinearModel with priors
         * \param X the X data
         * \param y the y data
         * \param rng a random-number generator class to use
         * \param init_beta the initial value of beta to use.  Defaults to 0 for all coefficients,
         * which, when combined with burn-in iterations, is probably fine.  If setting this affects
         * results, your burn-in period is probably too short.
         * \param init_s_sq the initial s^2 value to use.  Defaults to 1.0, which, when combined
         * with burn-in iterations, is probably fine.  If setting this affects results, your burn-in
         * period is probably too short.
         */
        Gibbs(LinearModel<K> &&lm,
                MatrixNK X, VectorN y,
                const RNG &rng,
                VectorK init_beta = VectorK::Zero(), double init_s_sq = 1.0)
        : model_{std::move(lm)}, X_{std::move(X)}, y_{std::move(y)}, rng_{rng},
            beta_{std::move(init_beta)}, h_{1.0 / init_s_sq},
            // Pre-perform various constant expressions:
            XtX_{(X_.transpose() * X_).eval()},
            Xty_{(X_.transpose() * y_).eval()},
            prior_Vinv_{model_.priorV().inverse().eval()},
            prior_Vinv_beta_{(prior_Vinv_ * model_.priorBeta()).eval()}
        {}

        /** Accesses the LinearModel attached to this Gibbs sampler.
         */
        LinearModel<K>& model() { return model_; }

        /** Performs `s` burn-in draws and discards them. */
        void burn(size_t s) {
            for (size_t i = 0; i < s; i++) {
                draw_one_(false);
            }
        }

        /** Performs `s` draws and stores them in the stored LinearModel, first clearing any old
         * draws.  It is highly recommended to call burn() before this method.
         */
        void draw(size_t S) {
            model_.draw_beta_ = Matrix<double, K, Dynamic>(K, S);
            model_.draw_s_sq_ = VectorXd(S);
            model_draw_ = 0;
            for (size_t s = 0; s < S; s++) {
                draw_one_(true);
            }
        }

    protected:
        /** Performs a draw of beta and then s^2 using the last (or initial) values of each, and
         * stores them as the most recent draw.  If the given parameter is true, they are also
         * stored in the model.
         */
        void draw_one_(bool store) {
            MatrixKK vbar = (prior_Vinv_ + h_ * XtX_).inverse();
            VectorK betabar = vbar * (prior_Vinv_beta_ + h_ * Xty_);
            MatrixKK Psi = vbar.llt().matrixL();
            VectorK z;
            std::normal_distribution<double> stdnormal;
            for (size_t i = 0; i < K; i++) {
                z[i] = stdnormal(rng_);
            }
            beta_ = betabar + Psi * z;

            VectorN resids = y_ - X_ * beta_;
            double s2bar = 1 / post_n_ * (resids.dot(resids) + model_.priorN() * model_.priorSSq());
            h_ = std::gamma_distribution<double>(post_n_/2.0, 2.0/(s2bar * post_n_))(rng_);

            if (store) {
                for (size_t k = 0; k < K; k++)
                    model_.draw_beta_(k, model_draw_) = beta_[k];
                model_.draw_s_sq_[model_draw_] = 1.0 / h_;
                model_draw_++;
            }
        }

    private:
        LinearModel<K> model_;
        MatrixNK X_;
        VectorN y_;
        RNG rng_;
        VectorK beta_;
        double h_;

        size_t model_draw_ = 0;

        /// The data-specific \f$X^\top X\f$
        MatrixKK XtX_;
        /// The data-specific \f$X^\top y\f$
        VectorK Xty_;
        /// prior V inverse, from the LinearModel
        MatrixKK prior_Vinv_;
        /// prior_Vinv_ times prior_beta, from the LinearModel
        VectorK prior_Vinv_beta_;
        /// New n size (= old n + N)
        const double post_n_ = model_.priorN() + N;
};

}
