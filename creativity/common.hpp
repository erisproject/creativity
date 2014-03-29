#pragma once
#include <eris/SharedMember.hpp>
#include <eris/Good.hpp>

/// \file creativity/common.hpp Storage class for global objects such as the money good.

namespace creativity {

/// The money good.  Set when created.
extern eris::SharedMember<eris::Good::Continuous> MONEY;

}
