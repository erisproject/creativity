#include "creativity/Creativity.hpp"
#include "creativity/state/FileStorage.hpp"
#ifndef CREATIVITY_SKIP_PGSQL
#include "creativity/state/PsqlStorage.hpp"
#endif
#include "creativity/state/MemoryStorage.hpp"
#include "creativity/belief/Demand.hpp"
#include "creativity/belief/Quality.hpp"
#include "creativity/belief/Profit.hpp"
#include <eris/Random.hpp>
#include <eris/intraopt/Callback.hpp>
#include <Eigen/Core>
#include <cmath>
#include <algorithm>
#include <regex>

namespace creativity {

using namespace creativity::state;
using namespace eris;
using namespace Eigen;

CreativitySettings& Creativity::set() {
    if (setup_sim_) throw std::logic_error("Cannot change creativity settings after setup()");
    else if (setup_read_) throw std::logic_error("Cannot change creativity settings after loading state data");
    return const_cast<CreativitySettings&>(parameters);
}

void Creativity::fileWrite(const std::string &filename) {
    if (setup_sim_) throw std::logic_error("Cannot call Creativity::fileWrite() after setup()");
    else if (setup_read_) throw std::logic_error("Cannot call Creativity::fileWrite() after loading state data");
    storage().first = Storage::create<FileStorage>(set_, filename, FileStorage::MODE::OVERWRITE);
}

void Creativity::fileRead(const std::string &filename) {
    if (setup_sim_) throw std::logic_error("Cannot call Creativity::fileRead() after setup()");
    storage().first = Storage::create<FileStorage>(set_, filename, FileStorage::MODE::READONLY);
    setup_read_ = true;
}

#ifndef CREATIVITY_SKIP_PGSQL
void Creativity::pgsql(std::string url, bool load_only, bool new_only) {
    if (setup_sim_) throw std::logic_error("Cannot call Creativity::pgsql() after setup()");

    // Figure out whether a "creativity=xyz" setting was provided, and if so, remove it
    int32_t load_id = -1;
    std::smatch results;
    if (std::regex_search(url, results, std::regex("([?&])creativity=(\\d+)(?=&|$)"))) {
        load_id = std::stoi(results[2].str());
        if (load_id <= 0) load_id = -1;
        url = results.prefix().str() + results.suffix().str();
    }

    if (load_id == -1 and setup_read_) throw std::logic_error("Cannot call Creativity::pgsql() after loading state data");
    else if (load_id > 0 and new_only) throw std::logic_error("Invalid postgresql URL given: creativity= option cannot be specified here");
    else if (load_id == -1 and load_only) throw std::logic_error("Invalid postgresql URL given: creativity= option must be specified here");

    storage().first = Storage::create<PsqlStorage>(set_, url, load_id);
    setup_read_ = load_id > 0;
}
#else
void Creativity::pgsql(std::string, bool, bool) {
    throw std::runtime_error("Postgresql support disabled at build time!");
}
#endif

void Creativity::checkParameters() {
#define PROHIBIT(FIELD, BAD) \
    if (parameters.FIELD BAD) throw std::domain_error("Invalid Creativity setting: parameters." #FIELD " " #BAD " is invalid")
    PROHIBIT(dimensions, < 1);
    PROHIBIT(readers, < 1);
    PROHIBIT(boundary, <= 0);
    PROHIBIT(book_distance_sd, < 0);
    PROHIBIT(book_quality_sd, < 0);
    PROHIBIT(reader_step_sd, < 0);
    PROHIBIT(reader_creation_shape, >= 1);
    PROHIBIT(reader_creation_scale_min, < 0);
    PROHIBIT(reader_creation_scale_max, < parameters.reader_creation_scale_min);
    PROHIBIT(cost_fixed, < 0);
    PROHIBIT(cost_unit, < 0);
    PROHIBIT(cost_piracy, < 0);
    PROHIBIT(income, < 0);
    PROHIBIT(piracy_link_proportion, < 0);
    PROHIBIT(piracy_link_proportion, > 1);
    PROHIBIT(initial.prob_write, < 0);
    PROHIBIT(initial.prob_write, > 1);
    PROHIBIT(initial.p_min, < 0);
    PROHIBIT(initial.p_max, < parameters.initial.p_min);
    PROHIBIT(initial.q_max, < parameters.initial.q_min);
    PROHIBIT(initial.q_min, < 0);
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
        if (not st.first) st.first = Storage::create<MemoryStorage>(set_);
    }

    std::uniform_real_distribution<double> unif_pmb(-parameters.boundary, parameters.boundary);
    std::uniform_real_distribution<double> unif_cr_shape(parameters.reader_creation_scale_min, parameters.reader_creation_scale_max);

    sim = Simulation::create();

    money = sim->spawn<Good::Continuous>();

    auto &rng = eris::Random::rng();

    if (parameters.piracy_link_proportion < 0 or parameters.piracy_link_proportion > 1)
        throw std::logic_error("Creativity piracy_link_proportion parameter is invalid (not in [0,1])");

    unsigned long max_links = parameters.readers * (parameters.readers - 1) / 2;
    unsigned long num_links = std::lround(parameters.piracy_link_proportion * max_links);

    // Track ids of created readers, to build the set of all possible edges as we go
    std::vector<eris_id_t> created;
    std::vector<std::pair<eris_id_t, eris_id_t>> potential_edges;
    if (num_links > 0) {
        created.reserve(parameters.readers);
        potential_edges.reserve(max_links);
    }

    for (unsigned int i = 0; i < parameters.readers; i++) {
        auto r = sim->spawn<Reader>(shared_from_this(),
                Position{unif_pmb(rng), unif_pmb(rng)},
                // (Nearly) non-informative priors for the rest:
                parameters.cost_fixed, parameters.cost_unit, parameters.income
                );
        r->writer_book_sd = parameters.book_distance_sd;
        r->writer_quality_sd = parameters.book_quality_sd;
        r->creation_shape = parameters.reader_creation_shape;
        r->creation_scale = unif_cr_shape(rng);

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

    sim->spawn<intraopt::FinishCallback>([this] { new_books_.clear(); });

    setup_sim_ = true;
    storage().first->updateSettings();
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

bool Creativity::sharing() const {
    if (!setup_sim_) throw std::logic_error("Cannot call sharing() on a non-live or unconfigured simulation");
    return parameters.piracy_begins > 0 and sim->t() >= parameters.piracy_begins;
}

double Creativity::priorWeight() const {
    if (!setup_sim_) throw std::logic_error("Cannot call priorWeight() on a non-live or unconfigured simulation");
    return sim->t() == parameters.piracy_begins ? parameters.prior_weight_piracy : parameters.prior_weight;
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
