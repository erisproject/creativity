#include "creativity/Policy.hpp"
#include <sstream>

namespace creativity {

Policy::Policy(const std::string &p) : Policy() {
    std::istringstream iss(p);
    std::string policy;
    while (std::getline(iss, policy, ',')) {
        if (policy == "public" || policy == "public-sharing")
            *this += PublicSharing();
        else if (policy == "vote" || policy == "voting" || policy == "public-voting")
            *this += PublicVoting();
        else if (policy == "catch" || policy == "catch-pirates")
            *this += CatchPirates();
        else if (policy == "none" || policy == "")
            /* ignore */;
        else
            throw std::runtime_error("Invalid/unknown Policy string: " + p);
    }
}

Policy& Policy::operator+=(const Policy &add) {
    if (add.publicSharing()) code_ |= POLICY_PUBLIC_SHARING;
    if (add.publicVoting())  code_ |= POLICY_PUBLIC_VOTING;
    if (add.catchPirates())  code_ |= POLICY_CATCH_PIRATES;
    return *this;
}

Policy& Policy::operator-=(const Policy &remove) {
    if (remove.publicSharing()) code_ &= ~POLICY_PUBLIC_SHARING;
    if (remove.publicVoting())  code_ &= ~POLICY_PUBLIC_VOTING;
    if (remove.catchPirates())  code_ &= ~POLICY_CATCH_PIRATES;
    return *this;
}



}
