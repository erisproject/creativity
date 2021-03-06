#include "creativity/CopyrightPolice.hpp"
#include "creativity/Reader.hpp"
#include <eris/Bundle.hpp>
#include <eris/random/util.hpp>

using namespace eris;
using namespace eris::random;

namespace creativity {

CopyrightPolice::CopyrightPolice(const Creativity &creativity) : creativity_{creativity} {
    if (tax() < 0) throw std::domain_error("CopyrightPolice creation error: lump sum tax cannot be negative");
    double mu = Creativity::evalPolynomial(tax(), creativity_.parameters.policy_catch_mu),
           sigma = Creativity::evalPolynomial(tax(), creativity_.parameters.policy_catch_sigma);
    if (not std::isfinite(mu)) throw std::domain_error("CopyrightPolice creation error: mu polynomial evaluated to a non-finite value");
    if (not std::isfinite(sigma) or sigma < 0) throw std::domain_error("CopyrightPolice creation error: sigma polynomial evaluated to a negative or non-finite value");

    normal_ = boost::math::normal(mu, sigma);
}

double CopyrightPolice::fine(unsigned pirated) {
    if (pirated == 0) return 0.0;

    double fine = Creativity::evalPolynomial(pirated, creativity_.parameters.policy_catch_fine);

    return fine > 0 ? fine : 0.0;
}

void CopyrightPolice::interApply() {
    auto lock = writeLock();
    if (tax() > 0) {
        Bundle tax_bill(creativity_.money, tax());
        for (auto &r : simulation()->agents<Reader>()) {
            lock.add(r);
            r->assets.transferApprox(tax_bill, 1e-6);
            lock.remove(r);
        }
    }
}

void CopyrightPolice::intraFinish() {
    // FIXME: record the number of people accused, and number guilty for analysis.  (Will need to
    // store them here, then extract them into the state object/file).
    //
    auto lock = writeLock();
    for (auto &r : simulation()->agents<Reader>()) {
        lock.add(r);
        // TODO: we don't currently remove readers from the simulation; if that ever changes, this
        // will need to be updated to do something with the fine for pirated books without living
        // authors.
        unsigned pirated = 0;
        std::unordered_map<SharedMember<Reader>, unsigned> author_shares;
        for (const auto &bk : r->newBooks()) {
            if (bk.second.get().pirated() and bk.first->author()->id() != id()) {
                ++pirated;

                ++author_shares[bk.first->author()];
            }
        }
        bool caught = rcoin(prob(pirated));
        bool guilty = pirated > 0;
        if (caught) {
            // FIXME: record this status (so that it is tracked in the reader state)
            r->alterUtility(-creativity_.parameters.policy_catch_cost);

            if (guilty) {
                double fine = Creativity::evalPolynomial(pirated, creativity_.parameters.policy_catch_fine);
                if (fine > 0) {
                    r->alterUtility(-fine);
                    for (const auto &a : author_shares) {
                        auto lock_a = lock.supplement(a.first);
                        a.first->alterUtility(a.second * fine / pirated);
                    }
                }
            }
        }
        lock.remove(r);
    }
}

}
