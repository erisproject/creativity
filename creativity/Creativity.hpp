#pragma once
#include <eris/Simulation.hpp>
#include <eris/Good.hpp>
#include <eris/Optimize.hpp>
#include <Eigen/Core>
#include "creativity/Book.hpp"
#include "creativity/CreativitySettings.hpp"

namespace creativity {

namespace state {
// Forward declaration
class Storage;
}

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

        /** Simulation parameters that are used to configure the simulation when calling setup().
         * To adjust these parameters, call set().
         */
        const CreativitySettings &parameters = set_;

        /** Provides writeable access to the parameters used to configure the simulation.  This
         * method may only be called before setup() or fileRead(); attempting to call it afterwards
         * will raise an exception.
         */
        CreativitySettings& set();

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

        /** Checks `.parameters` to make sure that all configured values are valid, throwing an
         * exception if any invalid values are found.
         *
         * This doesn't check that parameters are useful, it merely checks whether the parameters
         * can be used.  (For example, readers=1 is a useless setting, but is still permitted
         * because such a simulation can technically be created).
         *
         * This method is called automatically by setup(), but can be called otherwise as well.
         *
         * \throws std::domain_error if any of the values in `.parameters` are invalid.
         */
        void checkParameters();

        /** Creates the Simulation object and adds the configured number of readers plus some base
         * members (such as the common money good).  After calling this method, the simulation will
         * be available in `.sim`.
         *
         * \throws std::domain_error (via checkParameters()) if any parameters are invalid.
         */
        void setup();

        /** Static method that calculates a boundary given a number of readers, dimensions, and a
         * desired density.
         *
         * \throws std::logic_error if any of the parameters are <= 0.
         */
        static double boundaryFromDensity(uint32_t readers, uint32_t dimensions, double density);

        /** Static method that calculates a density given a number of readers, dimensions, and a
         * desired boundary location.
         *
         * \throws std::logic_error if any of the parameters are <= 0.
         */
        static double densityFromBoundary(uint32_t readers, uint32_t dimensions, double boundary);

        /** Returns the value of densityFromBoundary called with the current readers, dimensions,
         * and boundary settings.
         */
        double densityFromBoundary() const;

        /** Returns true if file sharing exists.  Attempt to call this on a Creativity object that
         * is not a live simulation, or is a live simulation but has not been set up yet, will raise
         * an exception.
         */
        bool sharing() const;

        /** Returns the prior multiplier currently in effect.  This is
         * CreativitySettings.prior_scale except in the first piracy period, during which it is
         * CreativitySettings.prior_scale_piracy, and in the initial burnin periods, during which it
         * is CreativitySettings.prior_scale_burnin.
         */
        double priorWeight() const;

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

        /** Stores the number of on-market books.  This is updated at the beginning of every new
         * simulation stage with the number of on-market books in the just-ended period.
         */
        unsigned long market_books = 0;

        /** Stores the lagged number of on-market books, that is, the number of books that were on
         * the market in the period before the period that just ended.  This is updated at the same
         * time as market_books.
         */
        unsigned long market_books_lagged = 0;

    protected:
        /* Default constructor is protected; construct by calling create(). */
        Creativity() = default;

    private:
        // True if this is a live simulation.  Exclusive of setup_read_.
        bool setup_sim_;

        // True if this has loaded a previously simulation data file.  Exclusive of setup_sim_.
        bool setup_read_;

        // The actual, non-const settings (the public `.parameters` is a const reference to this)
        CreativitySettings set_;

        /** Stores new books.  Read by Reader and added to by Book to avoid having the search the
         * entire set of Books for new books.  This is cleared at the end of every period.
         *
         * Obtain a lock by calling booksLock() before accessing.
         */
        std::vector<eris::SharedMember<Book>> new_books_;

        /** The storage object.  If not set by the time setup() is called, a MemoryStorage object
         * will be created and used.
         */
        std::shared_ptr<state::Storage> storage_;

        // Mutex controlling .storage_ and new_books_ access
        std::mutex storage_mutex_, new_books_mutex_;
};

}
