#!/bin/bash

set -e

for a in "$@"; do
    if [ "$a" == "--help" ]; then
        cat <<EOF
Usage: $0 [ARGUMENTS]

Invokes ./creativity-cli via ./creativity-random, using a suitable set of
randomized parameters.

Output is placed in ./results/DATE~~HASH/randomized/creativity-SEED.crstate.xz
where DATE and HASH are generated from the current git commit date and hash.
EOF
        exit 1;
    fi
done

dir="./results/$(git show -s --format=%cd~~%h @)/randomized"
mkdir -p "$dir"

exec ./creativity-random ./creativity-cli \
    --readers 'iU[100,200]' \
    --density 'U[0.25,4]' \
    --reader-step-mean 'U[.1,.9]' \
    --book-distance-mean '0' \
    --book-quality-sd '0' \
    --creation-fixed 'U[50,250]' \
    --creation-time 'iU[0,5]' \
    --reader-creation-scale-min 0 \
    --reader-creation-scale-range 'U[5,15]' \
    --cost-market 'U[0,50]' \
    --cost-unit 'U[0,10]' \
    --cost-piracy 'U[0,10]' \
    --initial-prob-write 'U[0.1,0.4]' \
    --initial-effort-min 'U[0,50]' \
    --initial-effort-range 'U[10,100]' \
    --initial-price-min 'U[0,10]' \
    --initial-price-range 'U[10,40]' \
    --initial-prob-keep 'U[0.25,0.75]' \
    --initial-keep-price 'U[0.25,0.75]' \
    --piracy-link-proportion 'U[0.05,0.25]' \
    --public-sharing-tax 'U[1,100]' \
    --threads 0 \
    --output "$dir/creativity-SEED.crstate" \
    "$@"
