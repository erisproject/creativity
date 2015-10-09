#!/bin/bash

# Invokes run-sim.sh many times.  Should be run from build dir.

if [ "$#" != 1 ] || ! [[ "$1" =~ ^[0-9]+$ ]]; then
    echo "Usage: $0 N -- invokes run-sim.sh N times"
    exit 1
fi

toexec="${0/many-sims.sh/run-sim.sh}"

for ((i = 1; i <= "$1"; i++)) do
    echo "Running simulation $i/$1"
    $toexec
done
