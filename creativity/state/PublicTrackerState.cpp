#include "creativity/state/PublicTrackerState.hpp"
#include "creativity/PublicTracker.hpp"

namespace creativity { namespace state {

PublicTrackerState::PublicTrackerState(const PublicTracker &pt)
    : id{pt.id()}, dl_tax{pt.dlTax()}, dl_unspent{pt.dlPool()}, vote_unspent{pt.votePool()}
{}

}}
