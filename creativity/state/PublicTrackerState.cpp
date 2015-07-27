#include "creativity/state/PublicTrackerState.hpp"

namespace creativity { namespace state {

PublicTrackerState::PublicTrackerState(const PublicTracker &pt)
    : id{pt.id()}, tax{pt.tax()}, unspent{pt.pool()}
{}

}}
