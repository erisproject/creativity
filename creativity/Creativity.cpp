#include "creativity/Creativity.hpp"
#include "creativity/state/FileStorage.hpp"
#include "creativity/PublicTracker.hpp"
#include "creativity/state/Storage.hpp"
#include "creativity/Reader.hpp"
#include "creativity/BookMarket.hpp"
#include <eris/Random.hpp>
#include <eris/Position.hpp>
#include <eris/interopt/Callback.hpp>
#include <eris/intraopt/Callback.hpp>
#include <cstddef>
#include <cmath>
#include <stdexcept>

namespace creativity {

using namespace creativity::state;
using namespace eris;

CreativitySettings& Creativity::set() {
    if (setup_sim_) throw std::logic_error("Cannot change creativity settings after setup()");
    else if (setup_read_) throw std::logic_error("Cannot change creativity settings after loading state data");
    return const_cast<CreativitySettings&>(parameters);
}

void Creativity::fileWrite(const std::string &filename) {
    write<FileStorage>(filename, FileStorage::MODE::OVERWRITE);
}

void Creativity::fileRead(const std::string &filename, bool xz_to_ram, bool copy_to_ram) {
    read<FileStorage>(filename, xz_to_ram, copy_to_ram);
}

void Creativity::checkParameters() {
#define PROHIBIT(FIELD, BAD) \
    if (parameters.FIELD BAD) throw std::domain_error("Invalid Creativity setting: parameters." #FIELD " " #BAD " is invalid")
    PROHIBIT(dimensions, < 1);
    PROHIBIT(readers, < 1);
    PROHIBIT(boundary, <= 0);
    PROHIBIT(book_distance_mean, < 0);
    PROHIBIT(book_quality_sd, < 0);
    PROHIBIT(reader_step_mean, < 0);
    PROHIBIT(reader_creation_shape, >= 1);
    PROHIBIT(reader_creation_scale_min, < 0);
    PROHIBIT(reader_creation_scale_range, < 0);
    PROHIBIT(cost_market, < 0);
    PROHIBIT(cost_unit, < 0);
    PROHIBIT(cost_piracy, < 0);
    PROHIBIT(income, < 0);
    PROHIBIT(piracy_link_proportion, < 0);
    PROHIBIT(piracy_link_proportion, > 1);
    PROHIBIT(initial.prob_write, < 0);
    PROHIBIT(initial.prob_write, > 1);
    PROHIBIT(initial.p_min, < 0);
    PROHIBIT(initial.p_range, < 0);
    PROHIBIT(initial.l_min, < 0);
    PROHIBIT(initial.l_range, < 0);
    PROHIBIT(initial.prob_keep, < 0);
    PROHIBIT(initial.prob_keep, > 1);
    PROHIBIT(initial.keep_price, < 0);
#undef PROHIBIT
}

void Creativity::setup() {
    if (setup_read_) throw std::logic_error("Cannot call Creativity::setup() after reading state data");
    else if (setup_sim_) throw std::logic_error("Creativity::setup() cannot be called twice");

    checkParameters();

    {
        auto st = storage();
        if (not st.first) st.first = Storage::create<FileStorage>(set_);
    }

    sim = Simulation::create();

    money = sim->spawn<Good::Continuous>();

    auto &rng = eris::Random::rng();
    std::uniform_real_distribution<double> unif_pmb(-parameters.boundary, parameters.boundary);
    std::uniform_real_distribution<double> unif_cr_shape(parameters.reader_creation_scale_min,
            parameters.reader_creation_scale_min + parameters.reader_creation_scale_range);

    Position initpos = Position::zero(parameters.dimensions);
    for (unsigned int i = 0; i < parameters.readers; i++) {
        // Draw a uniform value for each dimension
        for (size_t d = 0; d < parameters.dimensions; d++) initpos[d] = unif_pmb(rng);

        auto r = sim->spawn<Reader>(*this, initpos);

        r->creation_shape = parameters.reader_creation_shape;
        r->creation_scale = unif_cr_shape(rng);
    }


    sim->spawn<intraopt::FinishCallback>([this] { new_books_.clear(); });

    // Count the number of on-market books, store it, and update the lagged number of on-market
    // books to the previous number of on-market books.
    sim->spawn<interopt::BeginCallback>([this] {
        market_books_lagged = market_books;
        market_books = sim->countMarkets<BookMarket>();
    });

    setup_sim_ = true;
    storage().first->updateSettings();
}

void Creativity::run() {
    if (setup_read_) throw std::logic_error("Cannot call Creativity::run() on a preloaded simulation data file");
    else if (not setup_sim_) throw std::logic_error("Creativity::run() cannot be called before Creativity::setup()");


    // If piracy begins in the next period, set up the piracy network
    if (sim->t() + 1 == parameters.piracy_begins)
        createPiracyNetwork();

    // If the next period is the first when the PublicTracker becomes available, create it
    if (sim->t() + 1 == parameters.public_sharing_begins)
        sim->spawn<PublicTracker>(*this, parameters.public_sharing_tax);

    sim->run();

    storage().first->emplace_back(sim);
}

void Creativity::createPiracyNetwork() {

    if (parameters.piracy_link_proportion < 0 or parameters.piracy_link_proportion > 1)
        throw std::logic_error("Creativity piracy_link_proportion parameter is invalid (not in [0,1])");

    auto readers = sim->agents<Reader>();
    unsigned long max_links = readers.size() * (readers.size() - 1) / 2;
    unsigned long num_links = std::lround(parameters.piracy_link_proportion * max_links);

    // If piracy_link_proportion is 0, we don't actually have any piracy links
    if (num_links == 0) return;

    // Use reservoir sampling: the first num_links edges get added with probability 1.  Once filled
    // up, each new potential link makes the cut with probability num_links / i, where i is the
    // iteration number (beginning at 1).  If it makes the cut, it replaces a random element in the
    // current list of elements.
    //
    // If L is the number of links, the (L+1)th element has probability L/(L+1) of being added.  If
    // this happens, an arbitrary existing element has conditional probability 1/L of being removed
    // and so has unconditional probability of 1/(L+1) of being removed, which means an
    // unconditional probability of L/(L+1) of surviving (which, by design, equals the probability
    // of the new element being added).
    //
    // For i >= 2, the (L+i)th element has probability L/(L+i) of being incorporated.  If this
    // happens, an arbitrary element currently in the list has probability 1/L of being removed,
    // conditional on both being in the list in the first place, and on the new element being
    // incorporated--which are independent events.  Thus the probability of being replaced in this
    // round conditional on surviving the previous round is 1/L * L/(L+i) = 1/(L+i).  The
    // unconditional probability of being still in the list after element L+i is thus the
    // probability of being in the list in the previous step minus the probability of being replaced
    // by element L+i, or in other words, L/(L+i-1) - L/((L+i-1)(L+i)) = L/(L+i)--which again equals
    // the probability of the new element having been incorporated.
    //
    // Thus, by induction, after all N >= L steps, every element has probability L/N of being in the
    // list, and the list contains exactly L elements.

    auto &rng = eris::Random::rng();
    std::uniform_int_distribution<size_t> replace_dist(0, num_links-1);

    std::vector<std::pair<size_t, size_t>> edges;
    edges.reserve(num_links);
    unsigned long index = 1;
    for (size_t i = 1; i < readers.size(); i++) {
        for (size_t j = 0; j < i; j++, index++) {
            if (index <= num_links) {
                edges.emplace_back(i, j);
            }
            else {
                // Replace with probability num_links / index:
                if (std::uniform_int_distribution<size_t>(0, index)(rng) <= num_links) {
                    size_t replace = replace_dist(rng);
                    edges[replace].first = i;
                    edges[replace].second = j;
                }
            }
        }
    }

    for (auto &e : edges) {
        readers[e.first]->addFriend(readers[e.second]);
    }
}

double Creativity::boundaryFromDensity(uint32_t readers, uint32_t dimensions, double density) {
    if (readers == 0) throw std::logic_error("Cannot calculate boundary when readers == 0");
    if (dimensions == 0) throw std::logic_error("Cannot calculate boundary when dimensions == 0");
    if (density <= 0) throw std::logic_error("Cannot calculate boundary when density <= 0");

    const double r_over_d = readers/density;
    // Calculate the boundaries from the density.  Total hypervolume is (2*boundary)^(dimensions),
    // so to achieve `density` we need boundary set as the solution to:
    //     density = readers / ((2*boundary)^(dimensions))
    // which is:
    //     boundary = 1/2 * (readers / density)^(1/dimensions)
    // thus:
    return 0.5 *
        (dimensions == 1 ? r_over_d :
         dimensions == 2 ? std::sqrt(r_over_d) :
         dimensions == 3 ? std::cbrt(r_over_d) :
         std::pow(r_over_d, 1.0/dimensions));
}

double Creativity::densityFromBoundary(uint32_t readers, uint32_t dimensions, double boundary) {
    if (readers == 0) throw std::logic_error("Cannot calculate density when readers == 0");
    if (dimensions == 0) throw std::logic_error("Cannot calculate density when dimensions == 0");
    if (boundary <= 0) throw std::logic_error("Cannot calculate density when boundary <= 0");

    return readers / std::pow(2*boundary, dimensions);
}

double Creativity::densityFromBoundary() const {
    return densityFromBoundary(parameters.readers, parameters.dimensions, parameters.boundary);
}

bool Creativity::piracy() const {
    if (!setup_sim_) throw std::logic_error("Cannot call piracy() on a non-live or unconfigured simulation");
    return parameters.piracy_begins > 0 and sim->t() >= parameters.piracy_begins;
}

bool Creativity::publicSharing() const {
    if (!setup_sim_) throw std::logic_error("Cannot call publicSharing() on a non-live or unconfigured simulation");
    return parameters.public_sharing_begins > 0 and sim->t() >= parameters.public_sharing_begins;
}

double Creativity::priorWeight() const {
    if (!setup_sim_) throw std::logic_error("Cannot call priorWeight() on a non-live or unconfigured simulation");
    auto t = sim->t();
    return t == parameters.piracy_begins ? parameters.prior_scale_piracy :
        t == parameters.public_sharing_begins ? parameters.prior_scale_public_sharing :
        t <= parameters.burnin_periods ? parameters.prior_scale_burnin :
        parameters.prior_scale;
}

double Creativity::meanInitialQuality() const {
    const double &r = parameters.initial.l_range, &l = parameters.initial.l_min,
          &beta = parameters.reader_creation_shape;
    const double L = l+r;
    double val;
    if (beta == 0) {
        const double logLp1 = std::log(L+1), loglp1 = std::log(l+1);
        val = (L*logLp1 - L + logLp1)
            - (l*loglp1 - l + loglp1);
    }
    else {
        val = (std::pow(L+1, beta+1) / (beta+1) - L)
            - (std::pow(l+1, beta+1) / (beta+1) - l);
    }
    return val / r * (parameters.reader_creation_scale_min + parameters.reader_creation_scale_range / 2);
}

std::pair<std::vector<SharedMember<Book>>&, std::unique_lock<std::mutex>> Creativity::newBooks() {
    return std::pair<std::vector<SharedMember<Book>>&, std::unique_lock<std::mutex>>(
            new_books_, std::unique_lock<std::mutex>(new_books_mutex_));
}

std::pair<std::shared_ptr<Storage>&, std::unique_lock<std::mutex>> Creativity::storage() {
    return std::pair<std::shared_ptr<Storage>&, std::unique_lock<std::mutex>>(
            storage_, std::unique_lock<std::mutex>(storage_mutex_));
}


}
