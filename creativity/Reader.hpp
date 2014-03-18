#include <eris/Positional.hpp>
#include <eris/Agent.hpp>

namespace creativity {

class ReaderBase : public eris::Agent {
};

typedef eris::Positional<ReaderBase> Reader;

}
