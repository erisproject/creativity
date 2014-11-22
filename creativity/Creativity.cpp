#include "creativity/Creativity.hpp"
#include "creativity/state/FileStorage.hpp"
#include "creativity/state/MemoryStorage.hpp"
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

CreativitySettings& Creativity::set() {
    if (setup_sim_) throw std::logic_error("Cannot change creativity settings after setup()");
    else if (setup_read_) throw std::logic_error("Cannot change creativity settings after fileRead()");
    return const_cast<CreativitySettings&>(parameters);
}

void Creativity::fileWrite(const std::string &filename) {
    if (setup_sim_) throw std::logic_error("Cannot call Creativity::fileWrite() after setup()");
    else if (setup_read_) throw std::logic_error("Cannot call Creativity::fileWrite() after fileRead()");
    storage().first = std::make_shared<FileStorage>(filename, FileStorage::MODE::OVERWRITE, set_);
}

void Creativity::fileRead(const std::string &filename) {
    if (setup_sim_) throw std::logic_error("Cannot call Creativity::fileRead() after setup()");
    storage().first = std::make_shared<FileStorage>(filename, FileStorage::MODE::READONLY, set_);
    setup_read_ = true;
}

void Creativity::checkParameters() {
#define CHECK(FIELD, BAD) \
    if (parameters.FIELD BAD) throw std::domain_error("Invalid Creativity setting: parameters." #FIELD " " #BAD " is invalid")
    CHECK(dimensions, < 1);
    CHECK(readers, < 1);
    if (parameters.use_density) { CHECK(density, <= 0); }
    else { CHECK(boundary, <= 0); }
    CHECK(book_distance_sd, < 0);
    CHECK(book_quality_sd, < 0);
    CHECK(cost_fixed, < 0);
    CHECK(cost_unit, < 0);
    CHECK(cost_piracy, < 0);
    CHECK(income, < 0);
    CHECK(sharing_link_proportion, < 0);
    CHECK(sharing_link_proportion, > 1);
    CHECK(initial.prob_write, < 0);
    CHECK(initial.prob_write, > 1);
    CHECK(initial.p_min, < 0);
    CHECK(initial.p_max, < parameters.initial.p_min);
    // FIXME: Is q_min < 0 a problem?
    CHECK(initial.q_max, < parameters.initial.q_min);
    CHECK(initial.prob_keep, < 0);
    CHECK(initial.prob_keep, > 1);
    CHECK(initial.keep_price, < 0);
#undef CHECK
}

void Creativity::setup() {
    if (setup_read_) throw std::logic_error("Cannot call Creativity::setup() after reading a state file");
    else if (setup_sim_) throw std::logic_error("Creativity::setup() cannot be called twice");

    checkParameters();

    if (parameters.use_density) {
        set_.boundary = boundaryFromDensity(parameters.readers, parameters.dimensions, parameters.density);
        set_.use_density = false;
    }
    else {
        set_.density = densityFromBoundary(parameters.readers, parameters.dimensions, parameters.boundary);
    }

    storage().first = std::make_shared<MemoryStorage>(set_);


    std::uniform_real_distribution<double> unif_pmb{-parameters.boundary, parameters.boundary};

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
        r->writer_book_sd = parameters.book_distance_sd;
        r->writer_quality_sd = parameters.book_quality_sd;

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

bool Creativity::sharing() const {
    if (!setup_sim_) throw std::logic_error("Cannot call sharing() on a non-live or unconfigured simulation");
    return parameters.sharing_begins > 0 and sim->t() >= parameters.sharing_begins;
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
