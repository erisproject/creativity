#pragma once
#include <eris/Simulation.hpp>
#include <eris/Good.hpp>
#include <eris/Optimize.hpp>
#include <Eigen/Core>
#include "creativity/Book.hpp"
#include "creativity/state/MemoryStorage.hpp"

namespace creativity {

/** Central class for a creativity simulation; this class handles setting up the simulation
 * according to configured parameters and running the simulation.
 *
 * This class is the central code used by both the command line and GUI interfaces.
 */
class Creativity : private eris::noncopyable, public std::enable_shared_from_this<Creativity> {
    public:
        /** Creates a new Creativity object and returns a shared_ptr to it. */
        template<typename... Args>
        static std::shared_ptr<Creativity> create(Args&&... args) {
            return std::shared_ptr<Creativity>(new Creativity(std::forward<Args>(args)...));
        }

        /** The simulation object.  Will be emtpy until setup() is called. */
        std::shared_ptr<eris::Simulation> sim;

        /// The money good.  Will be empty until setup() is called.
        eris::SharedMember<eris::Good::Continuous> money;

        /** Simulation parameters that should be updated appropriately before calling setup().
         * These parameters may be read but should not be changed after setup() has been called.
         */
        struct {
            /// The number of readers in the simulation
            unsigned int readers = 100;

            /// The number of dimensions in the simulation
            unsigned int dimensions = 2;

            /** The density of the simulation in average number of readers per
             * \f$unit^{dimensions}\f$.  This, combined with `readers` and `dimensions` implicitly
             * defines the boundaries where the simulation space wraps to the opposite boundary.
             *
             * For example, for `readers = 100, dimensions = 2, density = 1`, the boundaries
             * will be at ±5 in each dimensions.  `readers = 100, dimensions = 3, density = 2` would
             * result in boundaries at ±1.842 in each of the three dimensions.
             */
            double density = 1.0;

            /** The standard deviation of a book.  An authored book will be located at a distance
             * drawn from \f$\left|N(0,s)\right|\f$, where \f$s\f$ is this value, in a random
             * direction from the author's location at the time of writing.
             *
             * The direction is drawn from a uniform distribution over the surface of the
             * hypersphere centred on the author with the randomly drawn radius.
             *
             * If this value is 0, the book is located exactly at the author's position at the time
             * of writing.
             */
            double book_distance_sd = 0.5;

            /** The standard deviation of a book quality draw.  When a reader obtains a book, his
             * subjective quality is drawn from \f$N(Q, s)\f$, where \f$Q\f$ is the book's base
             * quality as decided by the author and \f$s\f$ is this setting.
             *
             * As long as this setting is positive, this means book quality is subjective; if 0, all
             * readers perceive the book as having the same quality.
             */
            double book_quality_sd = 1.0;

            /** The fixed cost of keeping a book on the market.
             *
             * This cost can be configured per reader as well; this value only affects the initial
             * values when creating the simulation.
             *
             * \sa updateAllCosts(double, double)
             */
            double cost_fixed = 10.0;

            /** The unit cost of selling a copy of a book (note that only copies of books still on
             * the market can be sold).
             *
             * This cost can be configured per reader as well; this value only affects the initial
             * values when creating the simulation.
             *
             * \sa updateAllCosts(double, double)
             */
            double cost_unit = 1.0;

            /** The per-period external income readers receive.  Effort spent creating a book in a
             * period is subtracted from this amount.
             *
             * Like `cost_fixed` and `cost_unit` this only specifies the default reader income:
             * reader income can be updated on a individual reader level.
             */
            double income = 1000.0;

            /** The period in which the sharing network is introduced.
             */
            unsigned long sharing_begins = 100;

            /** The number of sharing/friendship links as a proportion of the maxinum number of
             * sharing links possible (which is \f$\frac{R(R-1)}{2}\f$, where \f$R\f$ is the number
             * of readers).
             *
             * Links as assigned randomly between agents when the simulation is initially set up,
             * with an equal probably of each possible link being selected.
             *
             * The default is 10% coverage of maximum potential links (rounded to the nearest
             * integer).  In the default 100-reader simulation, this is 495 links.
             *
             * The value must be in \f$[0, 1]\f$.
             */
            double sharing_link_proportion = 0.1;

        } parameters;

        /** Returns the simulation boundary.  If simulation setup is complete, this returns the
         * boundary calculated at the time the simulation was configured; otherwise it calculates
         * the boundary according to the current values of `parameters.readers`,
         * `parameters.dimensions`, and `parameters.density`.
         *
         * The value returned is the positive boundary coordinate, indicating boundaries at both
         * positive and negative values of the return value, in all dimensions.
         */
        double boundary() const;

        /** Store the creativity simulation results in the given file, overwriting any content the
         * file may already have.  The existing Storage object is released (and destroyed, unless
         * something else has copied its shared_ptr).  After calling this method, `.storage` will be
         * set to the new FileStorage object.
         *
         * This must be called before calling setup().
         *
         * The default simulation storage is an in-memory storage.
         */
        void fileWrite(const std::string &filename);

        /** Reads a creativity simulation storage file into storage().  The file is opened readonly,
         * and so any attempt to manipulate the simulation will fail.
         *
         * This method is typically called instead of calling setup() when loading a creativity
         * simulation record from disk.
         */
        void fileRead(const std::string &filename);

        /** Creates the Simulation object and adds the configured number of readers plus some base
         * members (such as the common money good).  After calling this method, the simulation will
         * be available in `.sim`.
         */
        void setup();

        /** Updates the costs of all simulation readers to the given values.  If setup() has not
         * been called, this updates cost_fixed and cost_unit.
         *
         * This method obtains a simulation run-lock, which means it cannot be invoked from within a
         * simulation run.
         *
         * \param cost_fixed the new fixed cost.  If negative, fixed costs are not changed.
         * \param cost_unit the new unit cost.  If negative, unit costs are not changed.
         */
        void updateAllCosts(double cost_fixed, double cost_unit = -1.0);

        /** Returns true if file sharing exists.  Attempt to call this on a Creativity object that
         * is not a live simulation, or is a live simulation but has not been set up yet, will raise
         * an exception.
         */
        bool sharing() const;

        /** Returns the period in which sharing becomes available in the current simulation.
         * Raises an exception if the Creativity object corresponds neither to a simulation state
         * nor live, setup simulation.
         */
        unsigned long sharingBegins() const;

        /** Establishes a lock on the new books storage and returns a pair consisting of the new
         * books reference and a unique lock on the books.  For optimal performance, store the
         * returned value in the smallest scope practical.
         *
         * \returns a pair with `.first` set to a reference to the new book vector and `.second` set
         * to the established lock.
         */
        std::pair<std::vector<eris::SharedMember<Book>>&, std::unique_lock<std::mutex>> newBooks();

        /** Establishes a lock on the storage object and returns a pair consisting of the storage
         * object and the unique lock on the object.  For optimal performance, store the returned
         * value in the smallest scope practical.
         *
         * \returns a pair with `.first` set to a reference to the storage shared_ptr and `.second`
         * set to the established lock.
         */
        std::pair<std::shared_ptr<state::Storage>&, std::unique_lock<std::mutex>> storage();

    protected:
        /* Default constructor is protected; construct by calling create(). */
        Creativity() = default;

    private:
        // True if this is a live simulation.  Exclusive of setup_read_.
        bool setup_sim_;

        // True if this has loaded a previously simulation data file.  Exclusive of setup_sim_.
        bool setup_read_;

        /** Stores new books.  Read by Reader and added to by Book to avoid having the search the
         * entire set of Books for new books.  This is cleared at the end of every period.
         *
         * Obtain a lock by calling booksLock() before accessing.
         */
        std::vector<eris::SharedMember<Book>> new_books_;

        /** The storage object.  Defaults to a MemoryStorage object, but can be replaced with a
         * FileStorage object by calling fileStorage().
         *
         * Obtain a lock by calling storageLock() before accessing.
         */
        std::shared_ptr<state::Storage> storage_ = std::make_shared<state::MemoryStorage>();

        // Mutex controlling .storage_ and new_books_ access
        std::mutex storage_mutex_, new_books_mutex_;

        // During setup() this is set to the wrapping boundary coordinates used by the simulation.
        double boundary_;

        // The period in which sharing first becomes available
        unsigned long sharing_begins_;

};

}
