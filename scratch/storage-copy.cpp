#include "creativity/state/FileStorage.hpp"
#include "creativity/state/Storage.hpp"
#include "creativity/CreativitySettings.hpp"
#include <eris/Random.hpp>
#include <eris/debug.hpp>

using namespace creativity;
using namespace creativity::state;
using namespace eris;
using namespace Eigen;

int main(int argc, char **argv) {
    if (argc != 3) {
        std::cerr << "Usage: storage-copy FROM TO -- reads FROM, copies it to TO\n";
        std::exit(1);
    }
    CreativitySettings cs;
    auto from = Storage::create<FileStorage>(cs, argv[1], FileStorage::MODE::READONLY);

    auto to = Storage::create<FileStorage>(cs, argv[2], FileStorage::MODE::OVERWRITE);

    for (auto &state : *from) {
        to->push_back(state);
    }
}
