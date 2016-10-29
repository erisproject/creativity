#pragma once
#include <eris/Good.hpp>
#include <eris/noncopyable.hpp>
#include "creativity/CreativitySettings.hpp"
#include "creativity/state/Storage.hpp"
#include <cstdint>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

namespace creativity { class Book; }

/// Primary namespace for all Creativity library code.
namespace creativity {

/** The bit for the public-sharing policy in `settings.policy`.
 * \internal
 */
constexpr uint32_t POLICY_PUBLIC_SHARING = 1 << 0;
/** The bit for the catch-pirates policy in `settings.policy`.
 * \internal
 */
constexpr uint32_t POLICY_CATCH_PIRATES  = 1 << 1;
/** The bit for the public-sharing-with-voting policy in `settings.policy`.
 * \internal
 */
constexpr uint32_t POLICY_PUBLIC_SHARING_VOTING = 1 << 2;

/** Central class for a creativity simulation; this class handles setting up the simulation
 * according to configured parameters and running the simulation.
 *
 * This class is the central code used by both the command line and GUI interfaces.
 */
class Creativity : private eris::noncopyable {
    public:
        /** Default constructor. */
        Creativity() = default;

        /** The simulation object.  Will be emtpy until setup() is called. */
        std::shared_ptr<eris::Simulation> sim;

        /// The money good.  Will be empty until setup() is called.
        eris::SharedMember<eris::Good> money;

        /** Simulation parameters that are used to configure the simulation when calling setup().
         * To adjust these parameters, call set().
         */
        const CreativitySettings &parameters{set_};

        /** Provides writeable access to the parameters used to configure the simulation.  This
         * method may only be called before setup() or fileRead(); attempting to call it afterwards
         * will raise an exception.
         */
        CreativitySettings& set();

        /** Creates a new Storage object with any type of StorageBackend class.  Writes any existing
         * and any new data to the new Storage object.
         *
         * Any existing simulation states are copied from the current Storage object, then the
         * existing Storage object is released (and destroyed, unless something else has copied its
         * shared_ptr).  After calling this method, `.storage` will be set to the new FileStorage
         * object.
         */
        template <class T, typename... Args>
        void write(Args&&... args) {
            auto new_storage = state::Storage::create<T>(set_, std::forward<Args>(args)...);
            auto spair = storage();
            auto &old_storage = spair.first;
            if (old_storage) {
                for (const auto &state : *old_storage) {
                    new_storage->push_back(state);
                    // Don't let the old storage get ahead of the new storage (we don't want to
                    // pointlessly slam a bunch of things into memory).
                    new_storage->flush(false);
                }
            }

            spair.first = new_storage;
        }

        /** Opens a creativity simulation storage source for reading via storage().  Typically
         * called instead of setup() when loading an existing creativity simulation record from a
         * stored location.
         *
         * This method cannot be combined with a call to setup(): while that method is used to
         * create a live simulation, this one reads a prior one.
         */
        template <class T, typename... Args>
        void read(Args&&... args) {
            if (setup_sim_) throw std::logic_error("Creativity: attempt to read() after calling setup()");
            storage().first = state::Storage::create<T>(set_, std::forward<Args>(args)...);
            setup_read_ = true;
        }

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

        /** Runs one iteration of the simulation.  This calls simulation()->run(), but also takes
         * care of other details (such as setting up the piracy network at the right time, creating
         * a public tracker agent when appropriate, and updating stored simulation state data).
         */
        void run();

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

        /** Returns true if piracy exists yet in the simulation.  Attempting to call this on a
         * Creativity object that is not a live simulation, or is a live simulation but has not been
         * set up yet, will raise an exception.
         */
        bool piracy() const;

        /** Returns true if a policy response exists yet in the simulation.  This must only be
         * called on a live, setup simulation.
         *
         * This returns false if the `policy_begins` period has not yet been reached, or if
         * there is no policy response configured at all.
         */
        bool policyActive() const;

        /** Returns true if public sharing is an enabled policy response *and* it is active in the
         * simulation (i.e. the simulation has reached the period in which the policy response
         * begins).  Attempting to call this on a Creativity object that is not a live simulation,
         * or is a live simulation but has not been set up yet, will raise an exception.
         */
        bool publicSharingActive() const;

        /** Returns true if public sharing with voting is an enabled policy response *and* it is
         * active in the simulation (i.e. the simulation has reached the period in which the policy
         * response begins).  Attempting to call this on a Creativity object that is not a live
         * simulation, or is a live simulation but has not been set up yet, will raise an exception.
         */
        bool publicVotingActive() const;

        /** Returns true if getting caught for piracy is a configured policy response *and* it is
         * active in the simulation (i.e. the simulation has reached the policy response stage).
         * Attempting to call this on a Creativity object that is not a live simulation, or is a
         * live simulation but has not been set up yet, will raise an exception.
         */
        bool catchPiratesActive() const;

        /** Returns the per-user policies taxes currently in effect.  Before the policy comes into
         * effect, this is 0; afterwards, it is whatever is the sum of the taxes of whichever
         * policies are enabled.
         *
         * \returns the quantity of the money good that must be surrended.
         */
        double policyTaxes() const;

        /** Returns the level of exogenous, disposable income, i.e. income minus policy taxes.  This
         * is simply a shortcut for `creativity.parameters.income - creativity.policyTaxes()`.  Note
         * that this does not include any income or costs associated with authorship.
         */
        double disposableIncome() const;

        /** Returns the prior multiplier currently in effect.  This is
         * CreativitySettings.prior_scale except in the first piracy period, during which it is
         * CreativitySettings.prior_scale_piracy, and in the initial burnin periods, during which it
         * is CreativitySettings.prior_scale_burnin.
         */
        double priorWeight() const;

        /** Returns the mean quality of books created with initial parameter distributions.
         *
         * With effort levels distributed \f$\ell \sim U[l,L]\f$, and creation_scale distributed
         * \f$\alpha \sim U[a,A]\f$, the expected value of the mean quality will be, since the two
         * distributions are independent, \f$E[\alpha]E[Q(\ell)] = \frac{A-a}{2} E[Q(\ell)]\f$,
         * where \f$Q(\ell)\f$ is the quality function not including the creation_scale multiplier.
         *
         * Integrating the quality function \f$q(\ell) = \alpha \frac{(\ell+1)^\beta - 1}{\beta}\f$
         * over the uniform range \f$[l,L]\f$ yields:
         * \f[
         *     E[Q(\ell)] = \frac{1}{L-l} \frac{1}{\beta} \left[\frac{(x+1)^{\beta+1}}{\beta+1} -
         *     x\right]_{x=l}^{L}
         * \f]
         * for \f$\beta \neq 0\f$ and
         * \f[
         *     E[Q(\ell)] = \frac{1}{L-l} \left[x \log(x+1) - x + \log(x+1)\right]_{x=l}^{L}
         * \f]
         * for the log version when \f$\beta = 0\f$.  Multiplying these expressions by
         * \f$\frac{1}{2}(A-a)\f$ gives the expected quality value of a book.
         *
         * Note that this function implicitly assumes that all readers' creation_scale and
         * creation_shape parameters remain set to the values as dictated by the current
         * `parameters` values.
         */
        double meanInitialQuality() const;

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

        /** Evaluates a polynomial at the value `x`.
         *
         * \param x the value at which to evaluate the polynomial.
         * \param it an iterator to the first coefficient of the polynomial.  The first value is
         * the constant, the second is the linear term, the \f$i\f$th term applies to the \f$x^i\f$
         * term.  The order of the polynomial is implicitly determined by the distance between
         * `begin` and `end`.  Typically `container.begin()`.
         * \param end the past-the-end iterator at which to stop evaluating the polynomial.  Typically
         * `container.end()`.
         */
        template <typename Iter, typename = typename std::enable_if<std::is_arithmetic<
            typename std::iterator_traits<Iter>::value_type>::value>::type>
        static double evalPolynomial(double x, Iter it, Iter end) {
            double v = 0;
            double xi = 1;
            for (; it != end; ++it) {
                v += *it * xi;
                xi *= x;
                if (xi == 0) break;
            }
            return v;
        }

        /** Evaluates a polynomial at the value `x`.  This is simply a shortcut for calling
         * `evalPolynomial(x, cont.begin(), cont.end())` for containers with an stl-like container
         * interface.
         */
        template <typename Container, typename = typename std::enable_if<std::is_arithmetic<
            typename Container::value_type>::value>::type>
        static double evalPolynomial(double x, const Container &cont) {
            return evalPolynomial(x, cont.begin(), cont.end());
        }

        /** Stores the number of on-market books.  This is updated at the beginning of every new
         * simulation stage (in inter-begin) with the number of on-market books in the just-ended period.
         */
        unsigned long market_books{0};

        /** Stores the average number of on-market books over the past `creation_time + 1` periods.
         * This is updated at the same time as market_books.
         */
        double market_books_avg{0.0};

    protected:
        /** Sets up the piracy network.  Called automatically by run() just before the
         * `piracy_begins` period.
         */
        void createPiracyNetwork();

    private:
        // True if this is a live simulation.  Exclusive of setup_read_.
        bool setup_sim_{false};

        // True if this has loaded a previously simulation data file.  Exclusive of setup_sim_.
        bool setup_read_{false};

        // The actual, non-const settings (the public `.parameters` is a const reference to this)
        CreativitySettings set_;

        /** Stores new books.  Read by Reader and added to by Book to avoid having the search the
         * entire set of Books for new books.  This is cleared at the end of every period.
         *
         * Obtain a lock by calling booksLock() before accessing.
         */
        std::vector<eris::SharedMember<Book>> new_books_;

        /** Linked list of number of market books over the last `creation_time + 1` periods, i.e. the
         * average number of books over a creativity cycle.
         */
        std::list<unsigned int> mkt_books_;

        /** The storage object.  If not set by the time setup() is called, a MemoryStorage object will
         * be created and used.
         */
        std::shared_ptr<state::Storage> storage_;

        // Mutex controlling .storage_ and new_books_ access
        std::mutex storage_mutex_, new_books_mutex_;
};

}
