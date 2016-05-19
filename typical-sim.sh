#!/bin/bash

set -e

toexec="./creativity-cli"
preargs=()
postargs=()
for a in "$@"; do
    if [ "$a" == "--help" ]; then
        cat <<EOF
Usage: $0 [ARGUMENTS]

Invokes ./creativity-cli (and so is typically run from the build directory)
using a default set of parameter values derived from the median values of
simulations that produce semi-consistent writing.  Any added ARGUMENTS (other
than --help and --dry-run) are added to the end of the ./creativity-cli
argument list, override the default.

Output is placed in

    ./results/\$DATE~~\$HASH/\$TYPICAL/creativity-\$SEED.crstate

where DATE and HASH are generated from the current git commit date and hash,
and TYPICAL defaults to "typical" if not externally set in the environment.

Example usage:

    TYPICAL=costly-piracy ./typical-sim.sh --piracy-cost 100

To see a list of the options that creativity-cli will be called with without
actually invoking it, add a --dry-run option.
EOF
        exit 1;
    elif [ "$a" == "--dry-run" ]; then
        preargs=("
Dry run; would have invoked:

$toexec")
        toexec="/bin/echo"
    else
        postargs+=("$a")
    fi
done

if [ -z "$TYPICAL" ]; then TYPICAL="typical"; fi
dir="./results/$(git show -s --date=short --format=%cd~~%h @)/${TYPICAL}"
mkdir -p "$dir"

exec $toexec "${preargs[@]}" \
    --readers 150 \
    --density 2.7 \
    --reader-step-mean 0.45 \
    --book-distance-mean 0 \
    --book-quality-sd 0 \
    --creation-fixed 150 \
    --creation-time 3 \
    --reader-creation-scale-min 0 \
    --reader-creation-scale-range 12.3 \
    --cost-market 26.4 \
    --cost-unit 5 \
    --cost-piracy 5 \
    --initial-prob-write 0.25 \
    --initial-effort-min 25 \
    --initial-effort-range 55 \
    --initial-price-min 5 \
    --initial-price-range 25 \
    --initial-prob-keep 0.5 \
    --initial-keep-price 0.5 \
    --piracy-link-proportion 0.15 \
    --public-sharing-tax 50 \
    --threads 0 \
    --output "$dir/creativity-SEED.crstate" \
    "${postargs[@]}"
