#pragma once
#include <Eigen/Core>

using Eigen::MatrixXd;
using Eigen::RowVectorXd;

namespace bayes {

/** Abstract base class for a PosteriorSimulator.  This class stores data matrices and provides stubs
 * for burning and generator posterior draws.
 *
 * Subclasses must provide, at a minimum, implementations of the draw() and defaultBurns() methods.
 */
class PosteriorSimulator {
    public:
        /// No default constructor
        PosteriorSimulator() = delete;

        /** Creates the base PosteriorSimulator with X and y data.
         *
         * \param X the X data to copy into the .X field.
         * \param y the y data to copy into the .y field.
         *
         * \throws std::runtime_error if X and y have different number of rows
         */
        PosteriorSimulator(const MatrixXd &X, const MatrixXd &y);

        /** Creates the base PosteriorSimulator with only X data.  y will be an empty matrix.
         *
         * \param X the X matrix to copy into the .X field.
         */
        PosteriorSimulator(const MatrixXd &X);

        /** Performs a draw from the posterior and returns it.  Subclasses must provide an
         * interpretation of the returned values (for example, for a normal linear model, the draw
         * might be model coefficients followed by an s^2 value).
         */
        virtual RowVectorXd draw() = 0;

        /** Performs `s` draws from the posterior and returns a matrix of the results.  Each row is
         * a draw as returned by draw().
         */
        virtual MatrixXd drawMany(unsigned int s);

        /** Performs a draw from the posterior and discards it.  This is used to burn away the
         * initial draws, and can also be used to take non-sequential draws from the posterior.
         *
         * The default implementation of this method simply calls draw() without storing the
         * returned result, but subclasses are free to override it to perform more efficient
         * versions where the result need not be built at all.
         */
        virtual void discard();

        /** Discards `b` draws from the posterior.  This is used to burn away the effect of initial
         * parameter values.  The default implementation simply called discard() the given number of
         * times.
         */
        virtual void burn(unsigned int b);

        /** The default number of burn draws for this simulator.  Different simulators have
         * different decent burn values, and must override this method with an appropriate amount.
         */
        virtual unsigned int defaultBurns() = 0;

        /** The "X" data for the model, generally n rows by k columns.  See the appropriate subclass
         * for the specific value expected in X.
         */
        MatrixXd X;

        /** The "y" data for the model, if appropriate.  Typically a column, but could be multiple
         * columns or empty depending on the needs of the underlying model.
         */
        MatrixXd y;
};

}
