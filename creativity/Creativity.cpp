#include "creativity/Creativity.hpp"
#include "creativity/state/FileStorage.hpp"
#include "creativity/belief/Demand.hpp"
#include "creativity/belief/Quality.hpp"
#include "creativity/belief/Profit.hpp"
#include <eris/Random.hpp>
#include <eris/intraopt/Callback.hpp>
#include <Eigen/Core>
#include <cmath>
#include <algorithm>

namespace creativity {

using namespace creativity::state;
using namespace eris;
using namespace Eigen;

double Creativity::boundary() const {
    if (setup_sim_ or setup_read_) return boundary_;

    if (parameters.readers == 0) throw std::logic_error("Cannot calculate boundary when parameters.readers == 0");
    if (parameters.dimensions == 0) throw std::logic_error("Cannot calculate boundary when parameters.dimensions == 0");
    if (parameters.density <= 0) throw std::logic_error("Cannot calculate boundary when density <= 0");

    // Calculate the boundaries from the density.  Total hypervolume is (2*boundary)^(dimensions),
    // so to achieve `density` we need boundary set as the solution to:
    //     density = readers / ((2*boundary)^(dimensions))
    // thus:
    return 0.5 *
        (parameters.dimensions == 1 ? parameters.readers/parameters.density :
         parameters.dimensions == 2 ? std::sqrt(parameters.readers/parameters.density) :
         parameters.dimensions == 3 ? std::cbrt(parameters.readers/parameters.density) :
         std::pow(parameters.readers/parameters.density, 1.0/parameters.dimensions));
}

void Creativity::fileWrite(const std::string &filename) {
    if (setup_sim_) throw std::logic_error("Cannot call Creativity::fileWrite() after setup()");
    else if (setup_read_) throw std::logic_error("Cannot call Creativity::fileWrite() after fileRead()");
    storage().first = std::make_shared<FileStorage>(filename, FileStorage::MODE::OVERWRITE);
}

void Creativity::fileRead(const std::string &filename) {
    if (setup_sim_) throw std::logic_error("Cannot call Creativity::fileRead() after setup()");
    auto st = storage();
    st.first = std::make_shared<FileStorage>(filename, FileStorage::MODE::READONLY);
    boundary_ = st.first->boundary();
    sharing_begins_ = st.first->sharingBegins();
    setup_read_ = true;
}

void Creativity::setup() {
    if (setup_read_) throw std::logic_error("Cannot call Creativity::setup() after reading a state file");
    else if (setup_sim_) throw std::logic_error("Creativity::setup() cannot be called twice");

    boundary_ = boundary();

    sharing_begins_ = parameters.sharing_begins;
    storage().first->sharingBegins(sharing_begins_);

    std::uniform_real_distribution<double> unif_pmb{-boundary_, boundary_};

    sim = Simulation::spawn();

    money = sim->create<Good::Continuous>();

    auto &rng = eris::Random::rng();

    if (parameters.sharing_link_proportion < 0 or parameters.sharing_link_proportion > 1)
        throw std::logic_error("Creativity sharing_link_proportion parameter is invalid (not in [0,1])");

    unsigned long max_links = parameters.readers * (parameters.readers - 1) / 2;
    unsigned long num_links = std::lround(parameters.sharing_link_proportion * max_links);

    // Track ids of created readers, to build the set of all possible edges as we go
    std::vector<eris_id_t> created;
    std::vector<std::pair<eris_id_t, eris_id_t>> potential_edges;
    if (num_links > 0) {
        created.reserve(parameters.readers);
        potential_edges.reserve(max_links);
    }

    for (unsigned int i = 0; i < parameters.readers; i++) {
        auto r = sim->create<Reader>(shared_from_this(),
                Position{unif_pmb(rng), unif_pmb(rng)},
                // (Nearly) non-informative priors for the rest:
                belief::Demand(parameters.dimensions, belief::Demand::parameters()),
                belief::Profit(parameters.dimensions, belief::Profit::parameters()),
                belief::Quality(belief::Quality::parameters()),
                parameters.cost_fixed, parameters.cost_unit, parameters.income
                );
        r->writer_book_sd = parameters.book_quality_sd;

        if (num_links > 0) {
            for (auto &other : created) {
                potential_edges.emplace_back(r->id(), other);
            }
            created.push_back(r->id());
        }
    }

    if (num_links > 0 and num_links < max_links) {
        // We aren't a full graph, so shuffle then throw away the ones we don't need
        auto &rng = Random::rng();
        std::shuffle(potential_edges.begin(), potential_edges.end(), rng);
        potential_edges.resize(num_links);
    }

    for (auto &e : potential_edges) {
        sim->agent<Reader>(e.first)->addFriend(sim->agent<Reader>(e.second));
    }

    sim->create<intraopt::FinishCallback>([this] { new_books_.clear(); });

    setup_sim_ = true;
}

void Creativity::updateAllCosts(double cost_fixed, double cost_unit) {
    if (setup_read_) throw std::logic_error("Cannot change costs: Creativity data is read-only");
    if (cost_fixed >= 0)
        parameters.cost_fixed = cost_fixed;
    if (cost_unit >= 0)
        parameters.cost_unit = cost_unit;
    else if (cost_fixed < 0)
        // Both are negative
        throw std::logic_error("updateAllCosts must have at least one non-negative cost to update");

    if (not setup_sim_) return;

    auto lock = sim->runLock();
    for (auto &r : sim->agents<Reader>()) {
        auto lock = r->writeLock();
        if (cost_fixed >= 0)
            r->cost_fixed = cost_fixed;
        if (cost_unit >= 0)
            r->cost_unit = cost_unit;
    }
}

bool Creativity::sharing() const {
    if (!setup_sim_) throw std::logic_error("Cannot call sharing() on a non-live or unconfigured simulation");
    return sim->t() >= sharing_begins_;
}

unsigned long Creativity::sharingBegins() const {
    if (not setup_sim_ and not setup_read_) throw std::logic_error("Cannot call sharingBegins() on an unconfigured Creativity object");
    return sharing_begins_;
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
