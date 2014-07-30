#include "creativity/Creativity.hpp"
#include "creativity/state/FileStorage.hpp"
#include <eris/Random.hpp>
#include <eris/intraopt/Callback.hpp>
#include <Eigen/Core>
#include <cmath>

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
    else if (setup_read_) throw std::logic_error("Cannot call Creativity::fileWriter() after fileRead()");
    storage().first = std::make_shared<FileStorage>(filename, FileStorage::MODE::OVERWRITE);
}

void Creativity::fileRead(const std::string &filename) {
    if (setup_sim_) throw std::logic_error("Cannot call Creativity::fileRead() after setup()");
    auto st = storage();
    st.first = std::make_shared<FileStorage>(filename, FileStorage::MODE::READONLY);
    boundary_ = st.first->boundary();
    setup_read_ = true;
    
}

void Creativity::setup() {
    if (setup_read_) throw std::logic_error("Cannot call Creativity::setup() after reading a state file");
    else if (setup_sim_) throw std::logic_error("Creativity::setup() cannot be called twice");

    boundary_ = boundary();

    std::uniform_real_distribution<double> unif_pmb{-boundary_, boundary_};

    sim = Simulation::spawn();

    money = sim->create<Good::Continuous>();

    auto &rng = eris::Random::rng();
    ERIS_DBG("Setting up readers");
    ERIS_DBG("FIXME: need to move belief parameters into .parameters");

    VectorXd demand_beta{8}; demand_beta << 0, -2, 0.5, -0.1, 0, 0, 0, 0;
    MatrixXd demand_V = MatrixXd::Identity(8, 8);
    double demand_s2 = 10, demand_n = 1;
    VectorXd profit_beta{5}; profit_beta << 0, 1, 0, 0, 0;
    MatrixXd profit_V = MatrixXd::Identity(5, 5);
    double profit_s2 = 10, profit_n = 1;
    VectorXd quality_beta{7}; quality_beta << 5, -1, 1, 0, 0, 0, 0.1;
    MatrixXd quality_V = MatrixXd::Identity(7, 7);
    double quality_s2 = 10, quality_n = 1;

    for (unsigned int i = 0; i < parameters.readers; i++) {
        auto r = sim->create<Reader>(shared_from_this(),
                Position{unif_pmb(rng), unif_pmb(rng)},
                belief::Demand{2, demand_beta, demand_s2, demand_V, demand_n},
                belief::Profit{2, profit_beta, profit_s2, profit_V, profit_n},
                belief::Quality{quality_beta, quality_s2, quality_V, quality_n},
                parameters.cost_fixed, parameters.cost_unit, parameters.income
                );
        r->writer_book_sd = parameters.book_quality_sd;
    }
    ERIS_DBG("Done with readers");

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

std::pair<std::vector<SharedMember<Book>>&, std::unique_lock<std::mutex>> Creativity::newBooks() {
    return std::pair<std::vector<SharedMember<Book>>&, std::unique_lock<std::mutex>>(
            new_books_, std::unique_lock<std::mutex>(new_books_mutex_));
}

std::pair<std::shared_ptr<Storage>&, std::unique_lock<std::mutex>> Creativity::storage() {
    return std::pair<std::shared_ptr<Storage>&, std::unique_lock<std::mutex>>(
            storage_, std::unique_lock<std::mutex>(storage_mutex_));
}


}
