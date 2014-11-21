#include "creativity/Reader.hpp"
#include "creativity/Creativity.hpp"
#include <eris/Simulation.hpp>

using namespace eris;
using namespace creativity;

int main() {
    auto creativity = Creativity::create();
    creativity->setup();

    auto r = creativity->sim->agents<Reader>().front();

    double q123 = r->creationQuality(12.3);
    double e123 = r->creationEffort(12.3);
    double qe123 = r->creationQuality(e123);
    double eq123 = r->creationEffort(q123);
    ERIS_DBGVAR(q123);
    ERIS_DBGVAR(e123);
    ERIS_DBGVAR(qe123);
    ERIS_DBGVAR(eq123);

    r->creation_scale = 15;
    q123 = r->creationQuality(12.3);
    e123 = r->creationEffort(12.3);
    qe123 = r->creationQuality(e123);
    eq123 = r->creationEffort(q123);
    ERIS_DBGVAR(q123);
    ERIS_DBGVAR(e123);
    ERIS_DBGVAR(qe123);
    ERIS_DBGVAR(eq123);

    r->creation_shape = 0;
    q123 = r->creationQuality(12.3);
    e123 = r->creationEffort(12.3);
    qe123 = r->creationQuality(e123);
    eq123 = r->creationEffort(q123);
    ERIS_DBGVAR(q123);
    ERIS_DBGVAR(e123);
    ERIS_DBGVAR(qe123);
    ERIS_DBGVAR(eq123);

    r->creation_shape = 0.5;
    q123 = r->creationQuality(12.3);
    e123 = r->creationEffort(12.3);
    qe123 = r->creationQuality(e123);
    eq123 = r->creationEffort(q123);
    ERIS_DBGVAR(q123);
    ERIS_DBGVAR(e123);
    ERIS_DBGVAR(qe123);
    ERIS_DBGVAR(eq123);
}
