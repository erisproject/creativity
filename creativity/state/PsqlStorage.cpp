#include "creativity/state/PsqlStorage.hpp"
#include "creativity/Creativity.hpp"
#include <eris/Random.hpp>
#include <cmath>
#include <sys/stat.h>
#include <boost/algorithm/string.hpp>

namespace creativity { namespace state {

using namespace creativity::belief;
using namespace Eigen;
using namespace eris;

PsqlStorage::PsqlStorage(const std::string &connect, unsigned int id) {

    conn_ = std::unique_ptr<pqxx::connection>(new pqxx::connection(connect));

    initialConnection();

    // Connected (on failure, the above throws)

    pqxx::work trans(*conn_);
    if (id == 0) {
        // New simulation: create a simulation row
        trans.exec("INSERT INTO simulation (seed) VALUES (" + std::to_string(eris::Random::seed()) + ")");
        id_ = trans.prepared("lastval").exec()[0][0].as<int32_t>();
    }
    else {
        auto result = trans.exec("SELECT * FROM simulation WHERE id = " + std::to_string(id));
        if (result.empty())
            throw std::out_of_range("Invalid id `" + std::to_string(id) + "' passed to PsqlStorage constructor: the given simulation is not in the database");
        id_ = id;
        seed_ = result[0]["seed"].as<int64_t>();
        have_settings = true;
    }
    trans.commit();
}

void PsqlStorage::initialConnection() {
    // Pg defaults to 15 digits of precision; increase this by 2, so that we get exact double values
    conn_->set_variable("extra_float_digits", "2");

    conn_->prepare("lastval", "SELECT lastval()");
    conn_->prepare("insert_setting", "INSERT INTO setting (sim, name, dbl, i64) VALUES ($1,$2,$3,$4)");
    conn_->prepare("update_setting", "UPDATE setting SET dbl = $3, i64 = $4 WHERE sim = $1 AND name = $2");
    conn_->prepare("get_settings", "SELECT name, dbl, i64 FROM setting WHERE sim = $1");
    conn_->prepare("num_states", "SELECT COUNT(*) FROM state WHERE sim = $1");
    conn_->prepare("insert_state", "INSERT INTO state (sim, t) VALUES ($1,$2)");
    conn_->prepare("select_state", "SELECT * FROM state WHERE sim = $1 AND t = $2");
    conn_->prepare("select_readers", "SELECT * FROM reader WHERE state = $1");
    conn_->prepare("select_books", "SELECT * FROM book WHERE state = $1");
    conn_->prepare("friend_ids", "SELECT friend_eris_id FROM friend WHERE reader = $1");
    conn_->prepare("select_library", "SELECT * FROM library WHERE reader = $1 ORDER BY book_eris_id");
    conn_->prepare("select_beliefs", "SELECT * FROM belief WHERE reader = $1");
    conn_->prepare("insert_reader", "INSERT INTO reader (state,eris_id,position,u,u_lifetime,cost_fixed,cost_unit,cost_piracy,income) VALUES "
                                                       "($1" ",$2"   ",$3"    ",$4,$5"     ",$6"      ",$7"     ",$8"       ",$9)");
    conn_->prepare("insert_friend", "INSERT INTO friend (reader,friend_eris_id) VALUES ($1,$2)");
    conn_->prepare("insert_library_book", "INSERT INTO library (reader,book_eris_id,type,new,quality) VALUES ($1,$2,$3,$4,$5)");
    conn_->prepare("insert_belief", "INSERT INTO belief (reader, type, k, noninformative, s2, n, beta, v_lower) VALUES ($1,$2,$3,$4,$5,$6,$7,$8)");
    conn_->prepare("insert_book", "INSERT INTO book (state,eris_id,author_eris_id,created,position,quality,price,revenue,revenue_lifetime,sales,sales_lifetime,pirated,pirated_lifetime,lifetime) VALUES "
                                                   "($1" ",$2"   ",$3"  ",$4"   ",$5"    ",$6"   ",$7" ",$8"   ",$9"            ",$10"",$11"         ",$12"  ",$13"           ",$14)");
}

void PsqlStorage::writeSettings(const CreativitySettings &settings) {
    std::unique_lock<std::mutex> lock(conn_mutex_);
    pqxx::work trans(*conn_);
    std::string prepst(have_settings ? "update_setting" : "insert_setting");
    trans.prepared(prepst)(id)("dimensions")()(settings.dimensions).exec();
    trans.prepared(prepst)(id)("readers")()(settings.readers).exec();
    trans.prepared(prepst)(id)("boundary")(settings.boundary)().exec();
    trans.prepared(prepst)(id)("book_distance_sd")(settings.book_distance_sd)().exec();
    trans.prepared(prepst)(id)("book_quality_sd")(settings.book_quality_sd)().exec();
    trans.prepared(prepst)(id)("reader_step_sd")(settings.reader_step_sd)().exec();
    trans.prepared(prepst)(id)("reader_creation_shape")(settings.reader_creation_shape)().exec();
    trans.prepared(prepst)(id)("reader_creation_scale_min")(settings.reader_creation_scale_min)().exec();
    trans.prepared(prepst)(id)("reader_creation_scale_max")(settings.reader_creation_scale_max)().exec();
    trans.prepared(prepst)(id)("cost_fixed")(settings.cost_fixed)().exec();
    trans.prepared(prepst)(id)("cost_unit")(settings.cost_unit)().exec();
    trans.prepared(prepst)(id)("cost_piracy")(settings.cost_piracy)().exec();
    trans.prepared(prepst)(id)("income")(settings.income)().exec();
    trans.prepared(prepst)(id)("piracy_begins")()(settings.piracy_begins).exec();
    trans.prepared(prepst)(id)("piracy_link_proportion")(settings.piracy_link_proportion)().exec();
    trans.prepared(prepst)(id)("prior_weight")(settings.prior_weight)().exec();
    trans.prepared(prepst)(id)("prior_weight_piracy")(settings.prior_weight_piracy)().exec();
    trans.prepared(prepst)(id)("initial.prob_write")(settings.initial.prob_write)().exec();
    trans.prepared(prepst)(id)("initial.q_min")(settings.initial.q_min)().exec();
    trans.prepared(prepst)(id)("initial.q_max")(settings.initial.q_max)().exec();
    trans.prepared(prepst)(id)("initial.p_min")(settings.initial.p_min)().exec();
    trans.prepared(prepst)(id)("initial.p_max")(settings.initial.p_max)().exec();
    trans.prepared(prepst)(id)("initial.prob_keep")(settings.initial.prob_keep)().exec();
    trans.prepared(prepst)(id)("initial.keep_price")(settings.initial.keep_price)().exec();
    trans.prepared(prepst)(id)("initial.belief_threshold")()(settings.initial.belief_threshold).exec();
    trans.commit();
    dimensions_ = settings.dimensions;
    have_settings = true;
}

void PsqlStorage::readSettings(CreativitySettings &settings) const {
    std::unique_lock<std::mutex> lock(conn_mutex_);
    if (not have_settings) return;
    pqxx::work trans(*conn_);
    auto res = trans.prepared("get_settings")(id).exec();
    int found = 0;
    for (auto row : res) {
        found++;
        std::string setting(row[0].c_str());
        if (setting == "dimensions") settings.dimensions = row[2].as<uint32_t>();
        else if (setting == "readers") settings.readers = row[2].as<uint32_t>();
        else if (setting == "boundary") settings.boundary = row[1].as<double>();
        else if (setting == "book_distance_sd") settings.book_distance_sd = row[1].as<double>();
        else if (setting == "book_quality_sd") settings.book_quality_sd = row[1].as<double>();
        else if (setting == "reader_step_sd") settings.reader_step_sd = row[1].as<double>();
        else if (setting == "reader_creation_shape") settings.reader_creation_shape = row[1].as<double>();
        else if (setting == "reader_creation_scale_min") settings.reader_creation_scale_min = row[1].as<double>();
        else if (setting == "reader_creation_scale_max") settings.reader_creation_scale_max = row[1].as<double>();
        else if (setting == "cost_fixed") settings.cost_fixed = row[1].as<double>();
        else if (setting == "cost_unit") settings.cost_unit = row[1].as<double>();
        else if (setting == "cost_piracy") settings.cost_piracy = row[1].as<double>();
        else if (setting == "income") settings.income = row[1].as<double>();
        else if (setting == "piracy_begins") settings.piracy_begins = row[2].as<eris_time_t>();
        else if (setting == "piracy_link_proportion") settings.piracy_link_proportion = row[1].as<double>();
        else if (setting == "prior_weight") settings.prior_weight = row[1].as<double>();
        else if (setting == "prior_weight_piracy") settings.prior_weight_piracy = row[1].as<double>();
        else if (setting == "initial.prob_write") settings.initial.prob_write = row[1].as<double>();
        else if (setting == "initial.q_min") settings.initial.q_min = row[1].as<double>();
        else if (setting == "initial.q_max") settings.initial.q_max = row[1].as<double>();
        else if (setting == "initial.p_min") settings.initial.p_min = row[1].as<double>();
        else if (setting == "initial.p_max") settings.initial.p_max = row[1].as<double>();
        else if (setting == "initial.prob_keep") settings.initial.prob_keep = row[1].as<double>();
        else if (setting == "initial.keep_price") settings.initial.keep_price = row[1].as<double>();
        else if (setting == "initial.belief_threshold") settings.initial.belief_threshold = row[2].as<int32_t>();
        else throw std::out_of_range("Simulation with id `" + std::to_string(id) + "' has unknown setting `" + setting + "'");
    }
    if (found != 25) throw std::out_of_range("Simulation with id `" + std::to_string(id) + "' has invalid settings");
    dimensions_ = settings.dimensions;
    trans.commit();
}

size_t PsqlStorage::size() const {
    std::unique_lock<std::mutex> lock(conn_mutex_);
    pqxx::work trans(*conn_);
    size_t s = trans.prepared("num_states")(id).exec()[0][0].as<size_t>();
    trans.commit();
    return s;
}

void PsqlStorage::thread_insert(std::shared_ptr<const State> &&state) {
    std::unique_lock<std::mutex> lock(conn_mutex_);
    pqxx::work trans(*conn_);
    trans.prepared("insert_state")(id)(state->t).exec();
    int32_t state_id = trans.prepared("lastval").exec()[0][0].as<int32_t>();

    for (auto &r : state->readers) {
        insertReader(r.second, state_id, trans);
    }
    for (auto &bp : state->books) {
        auto &b = bp.second;
        trans.prepared("insert_book")(state_id)(b.id)(b.author)(b.created)
            (createDoubleArray(b.position.begin(), b.position.end()))(b.quality)(b.price)
            (b.revenue)(b.revenue_lifetime)(b.sales)(b.sales_lifetime)
            (b.pirated)(b.pirated_lifetime)(b.lifetime).exec();
    }
    trans.commit();
}

std::shared_ptr<const State> PsqlStorage::load(eris_time_t t) const {
    std::unique_lock<std::mutex> lock(conn_mutex_);
    std::shared_ptr<const State> ret;
    State *st_ptr = new State();
    State &state = *st_ptr;
    ret.reset(st_ptr);

    pqxx::work trans(*conn_);
    auto state_row = trans.prepared("select_state")(id)(t).exec();
    if (state_row.empty()) throw std::out_of_range("state::PsqlStorage: requested State index is invalid");

    int state_id = state_row[0]["id"].as<int>();
    state.t = state_row[0]["t"].as<eris_time_t>();

    for (auto reader_row : trans.prepared("select_readers")(state_id).exec()) {
        state.readers.insert(readReader(reader_row, trans));
    }
    for (auto book_row : trans.prepared("select_books")(state_id).exec()) {
        state.books.insert(readBook(book_row, trans));
    }
    trans.commit();

    return ret;
}


void PsqlStorage::insertReader(const ReaderState &reader, int32_t state, pqxx::work &trans) {

    trans.prepared("insert_reader")(state)(reader.id)(createDoubleArray(reader.position.begin(), reader.position.end()))
            (reader.u)(reader.u_lifetime)(reader.cost_fixed)(reader.cost_unit)(reader.cost_piracy)(reader.income).exec();

    int64_t rowid = trans.prepared("lastval").exec()[0][0].as<int64_t>();

    // Friends:
    for (auto &f : reader.friends) {
        trans.prepared("insert_friend")(rowid)(f).exec();
    }

    // Library:
    for (auto &l : reader.library) {
        auto &bid = l.first;
        auto insert = trans.prepared("insert_library_book");
        insert(rowid)(bid);
        if (reader.wrote.count(bid) > 0) insert("wrote")(false);
        else insert(reader.library_purchased.count(bid) > 0 ? "bought" : "pirated")(reader.new_books.count(bid) > 0);
        insert(l.second);
        insert.exec();
    }

    // Beliefs:
    insertBelief(rowid, "profit", reader.profit, trans);
    insertBelief(rowid, "profit_extrap", reader.profit_extrap, trans);
    insertBelief(rowid, "demand", reader.demand, trans);
    insertBelief(rowid, "quality", reader.quality, trans);
    for (auto &b : reader.profit_stream) {
        if (b.second.K() > 0) insertBelief(rowid, "profit_stream", b.second, trans);
    }
}


void PsqlStorage::insertBelief(eris_id_t dbid, const std::string &type, const Linear &belief, pqxx::work &trans) {
    auto insert = trans.prepared("insert_belief")(dbid)(type)(belief.K())(belief.noninformative());
    if (belief.noninformative())
        insert()()()();
    else {
        insert(belief.s2())(belief.n());
        std::vector<double> coef;
        auto &K = belief.K();
        coef.reserve(K*(K+1)/2);
        coef.resize(belief.K());
        for (unsigned int k = 0; k < K; k++) {
            coef[k] = belief.beta()[k];
        }
        insert(createDoubleArray(coef.cbegin(), coef.cend()));
        coef.resize(K*(K+1)/2);
        for (unsigned int r = 0, i = 0; r < K; r++) for (unsigned int c = 0; c <= r; c++, i++)
            coef[i] = belief.V()(r,c);
        insert(createDoubleArray(coef.cbegin(), coef.cend()));
    }
    insert.exec();
}

std::pair<eris::eris_id_t, ReaderState> PsqlStorage::readReader(const pqxx::tuple &reader_row, pqxx::work &trans) const {

    int64_t id = reader_row["id"].as<int64_t>(); // The internal id (NOT the eris_id_t)

    auto pair = std::make_pair<eris_id_t, ReaderState>(
            reader_row["eris_id"].as<eris_id_t>(),
            ReaderState(dimensions_));

    ReaderState &r = pair.second;
    r.id = pair.first;

    r.position = parseDoubleArray(reader_row["position"].c_str());
#define READER_READDBL(P) r.P = reader_row[#P].as<double>()
    READER_READDBL(u);
    READER_READDBL(u_lifetime);
    READER_READDBL(cost_fixed);
    READER_READDBL(cost_unit);
    READER_READDBL(cost_piracy);
    READER_READDBL(income);
#undef READER_READDBL

    // Friends:
    for (auto f : trans.prepared("friend_ids")(id).exec()) {
        r.friends.insert(f[0].as<eris_id_t>());
    }

    // Library:
    for (auto l : trans.prepared("select_library")(id).exec()) {
        eris_id_t book_id = l["book_eris_id"].as<eris_id_t>();
        r.library.emplace(book_id, l["quality"].as<double>()); // Add to library
        bool new_book = l["new"].as<bool>();

        if (new_book) r.new_books.insert(book_id); // Add to new_books

        std::string type(l["type"].c_str());
        if (type == "wrote")
            r.wrote.insert(book_id); // Add to wrote
        else if (type == "bought") {
            r.library_purchased.insert(book_id); // Add to library_purchased
            if (new_book) r.new_purchased.insert(book_id); // Add to new_purchased
        }
        else if (type == "pirated") {
            r.library_pirated.insert(book_id); // Add to library_pirated
            if (new_book) r.new_pirated.insert(book_id); // Add to new_pirated
        }
        else throw std::out_of_range("Reader with database id `" + std::to_string(id) + "' has library row with an invalid type");
    }

    // Beliefs:
    for (auto b : trans.prepared("select_beliefs")(id).exec()) {
        int k = b["k"].as<int>();
        if (k <= 0) continue; // Shouldn't happen, but if it does, the default-constructed beliefs are fine
        bool noninformative = b["noninformative"].as<bool>();
        std::string type(b["type"].c_str());

        if (noninformative) {
            // If noninformative, s2, n, beta, etc. are null
            if (type == "profit") r.profit = Profit(dimensions_, k);
            else if (type == "profit_extrap") r.profit_extrap = Profit(dimensions_, k);
            else if (type == "demand") r.demand = Demand(dimensions_, k);
            else if (type == "quality") r.quality = Quality(k);
            else if (type == "profit_stream") r.profit_stream.emplace((unsigned int) k, ProfitStream(k));
            else throw std::out_of_range("Reader with database id `" + std::to_string(id) + "' has invalid belief type `" + type + "'");
        }
        else {
            double s2 = b["s2"].as<double>();
            double n = b["n"].as<double>();
            auto beta = parseDoubleArray(b["beta"].c_str(), k);
            auto V = parseDoubleArray(b["v_lower"].c_str(), k*(k+1)/2);
            if (type == "profit") r.profit = Profit(dimensions_, beta, s2, V, n);
            else if (type == "profit_extrap") r.profit_extrap = Profit(dimensions_, beta, s2, V, n);
            else if (type == "demand") r.demand = Demand(dimensions_, beta, s2, V, n);
            else if (type == "quality") r.quality = Quality(beta, s2, V, n);
            else if (type == "profit_stream") r.profit_stream.emplace((unsigned int) k, ProfitStream(beta, s2, V, n));
            else throw std::out_of_range("Reader with database id `" + std::to_string(id) + "' has invalid belief type `" + type + "'");
        }
    }

    return pair;
}

std::pair<eris::eris_id_t, BookState> PsqlStorage::readBook(const pqxx::tuple &book_row, pqxx::work &) const {
    auto pair = std::make_pair<eris_id_t, BookState>(
            book_row["eris_id"].as<eris_id_t>(),
            BookState(dimensions_));

    BookState &b = pair.second;
    b.id = pair.first;
    b.position = parseDoubleArray(book_row["position"].c_str());
    auto price = book_row["price"];
    b.price = price.is_null() ? std::numeric_limits<double>::quiet_NaN() : price.as<double>();
    b.author = book_row["author_eris_id"].as<double>();
#define BOOK_READ(P, T) b.P = book_row[#P].as<T>()
    BOOK_READ(quality, double);
    BOOK_READ(revenue, double);
    BOOK_READ(revenue_lifetime, double);
    BOOK_READ(sales, unsigned int);
    BOOK_READ(sales_lifetime, unsigned int);
    BOOK_READ(pirated, unsigned int);
    BOOK_READ(pirated_lifetime, unsigned int);
    BOOK_READ(lifetime, unsigned int);
    BOOK_READ(created, eris_time_t);
#undef BOOK_READ

    return pair;
}


std::vector<double> PsqlStorage::parseDoubleArray(const std::string &doubles, int expected_size) {
    if (doubles.front() != '{' or doubles.back() != '}')
        throw std::invalid_argument("Invalid argument to parseDoubleArray: array delimiters not found in string");

    std::string csv(doubles.substr(1, doubles.length()-2));

    std::vector<std::string> dbl_strs;

    boost::split(dbl_strs, csv, [](const char &c) -> bool { return c == ','; });
    if (expected_size >= 0 and dbl_strs.size() != (unsigned) expected_size)
        throw std::invalid_argument("Invalid argument to parseDoubleArray: array has incorrect size");

    std::vector<double> dbls;
    dbls.resize(dbl_strs.size());
    int i = 0;
    for (auto &s : dbl_strs) pqxx::from_string(s.c_str(), dbls[i++]); // Might throw if invalid data found
    return dbls;
}

std::pair<pqxx::connection&, std::unique_lock<std::mutex>> PsqlStorage::connection() {
    return std::pair<pqxx::connection&, std::unique_lock<std::mutex>>(
            std::ref(*conn_),
            std::unique_lock<std::mutex>(conn_mutex_));
}

}}
