#include "creativity/belief/Linear.hpp"
#include "creativity/state/State.hpp"
#include "creativity/state/FileStorage.hpp"
#include <eris/Random.hpp>
#include <eris/debug.hpp>

using namespace creativity::belief;
using namespace creativity::state;
using namespace eris;
using namespace Eigen;

// Test class to exposed protected fields of FileStorage
class FS : public FileStorage {
    public:
        FS(std::string filename) : FileStorage(filename, MODE::OVERWRITE) {}

        std::fstream &f = f_;

        void test_write_model(Linear &l) {
            writeBelief(l);
        }

        belief_data test_read_model() {
            return readBelief();
        }

        unsigned int header_size() { return HEADER::size; }
};

auto rng = Random::rng();
std::uniform_real_distribution<double> r_pm10(-10, 10), r_01(0, 1);

Linear random_linear(int k) {
    VectorXd beta(k);
    for (int i = 0; i < k; i++) beta[i] = r_pm10(rng);
    MatrixXd V(k,k);
    for (int r = 0; r < k; r++) {
        for (int c = r; c < k; c++) {
            double v = r_01(rng);
            V(r,c) = v;
            V(c,r) = v;
        }
    }

    return Linear(beta, 100*r_01(rng), V, 10+100*r_01(rng));
}

int main() {
    Linear lin_default1, lin_default2;
    Linear lin5_noninf(5), lin1_noninf(1), lin55_noninf(55);
    Linear lin5_rand = random_linear(5),
           lin2_rand = random_linear(2),
           lin95_rand = random_linear(95);
    Linear lin_boring(VectorXd::Zero(7), 1.0, MatrixXd::Identity(7, 7), 25);

    // Write out a test file with each of the models to allow easy verification of a file containing
    // a single model:
    FS("storage.test.dc1").test_write_model(lin_default1);
    FS("storage.test.dc2").test_write_model(lin_default2);
    FS("storage.test.ni1").test_write_model(lin1_noninf);
    FS("storage.test.ni5").test_write_model(lin5_noninf);
    FS("storage.test.ni55").test_write_model(lin55_noninf);
    FS("storage.test.r2").test_write_model(lin2_rand);
    FS fsr5("storage.test.r5");
    fsr5.test_write_model(lin5_rand);
    fsr5.test_write_model(lin2_rand);
    fsr5.test_write_model(lin5_rand);
    fsr5.test_write_model(lin95_rand);
    fsr5.test_write_model(lin95_rand);
    fsr5.test_write_model(lin95_rand);
    fsr5.test_write_model(lin95_rand);
    fsr5.test_write_model(lin95_rand);
    fsr5.test_write_model(lin95_rand);
    fsr5.test_write_model(lin2_rand);
    fsr5.test_write_model(lin5_rand);
    FS("storage.test.r95").test_write_model(lin95_rand);

#define PRINTMODEL(m) std::cout << #m ": " << m << "\n    "; if (m.K()==0) std::cout << "(default constructed)"; else { if (m.noninformative()) std::cout << "(noninformative), "; std::cout << "n=" << m.n() << ", s2=" << m.s2() << ", V=\n" << m.V() << "\n"; }
    PRINTMODEL(lin_default1);
    PRINTMODEL(lin_default2);
    PRINTMODEL(lin1_noninf);
    PRINTMODEL(lin5_noninf);
//    PRINTMODEL(lin55_noninf);
    PRINTMODEL(lin2_rand);
    PRINTMODEL(lin5_rand);
//    PRINTMODEL(lin95_rand);

    bool test_pass = true;

#define TESTEQ(V1, V2, M) ERIS_DBGVAR(V1.M == V2.M()); if (test_pass and not (V1.M == V2.M())) test_pass = false;
#define COMPAREMODEL(N, LOC) { \
    fsr5.f.seekg(LOC); \
    auto read##N = fsr5.test_read_model(); \
    ERIS_DBG("Checking model at " #LOC); \
    TESTEQ(read##N, lin##N##_rand, beta); \
    TESTEQ(read##N, lin##N##_rand, n); \
    TESTEQ(read##N, lin##N##_rand, s2); \
    TESTEQ(read##N, lin##N##_rand, V); }

    COMPAREMODEL(5, 0x200);
    COMPAREMODEL(2, 0x2bc);
    COMPAREMODEL(95, 0x308);

    COMPAREMODEL(5, 0x300);
    COMPAREMODEL(5, 0x94cc);

    COMPAREMODEL(2, 0x94c4);

    COMPAREMODEL(95, 0x949c);
    COMPAREMODEL(95, 0x94a4);
    COMPAREMODEL(95, 0x94ac);
    COMPAREMODEL(95, 0x94b4);
    COMPAREMODEL(95, 0x94bc);

    if (test_pass) std::cout << "Tests passed\n\n";
    else std::cout << "Tests failed!  (compile with -DERIS_DEBUG if no debugging output above)\n\n";
}
