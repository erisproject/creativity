#include "creativity/state/PublicTrackerState.hpp"
#include "creativity/PublicTracker.hpp"

namespace creativity { namespace state {

PublicTrackerState::PublicTrackerState(const PublicTracker &pt)
    : id{pt.id()}, tax{pt.tax()}, unspent{pt.pool()}
{}

}}
