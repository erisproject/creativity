#!/bin/bash

set -e

for a in "$@"; do
    if [ "$a" == "--help" ]; then
        cat <<EOF
Usage: $0 [ARGUMENTS]

Invokes ./creativity-cli using a fixed set of parameter values derived from
the median values of simulations that produce consistent writing.

Output is placed in
    ./results/\${TYPICAL}~~\${DATE}~~\${HASH}/creativity-\${SEED}.crstate
where DATE and HASH are generated from the current git commit date and hash,
and TYPICAL defaults to "typical" if not externally set.
EOF
        exit 1;
    fi
done

if [ -z "$TYPICAL" ]; then TYPICAL="typical"; fi
dir="./results/${TYPICAL}~~$(git show -s --format=%cI~~%h @)"
mkdir -p "$dir"

exec ./creativity-cli \
    --readers 150 \
    --density 2.8 \
    --reader-step-mean 0.45 \
    --book-distance-mean 0 \
    --book-quality-sd 0 \
    --creation-fixed 150 \
    --creation-time 4 \
    --reader-creation-scale-min 0 \
    --reader-creation-scale-range 12.6 \
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
    "$@"
