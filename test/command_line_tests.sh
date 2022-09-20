#!/usr/bin/env bash

set -euo pipefail

julia --project=. deps/build.jl
OUT=$?
if [ "$OUT" != 0 ]; then
    echo build fail
    exit 1
else
    echo build OK
fi 

export PATH="$HOME/.julia/bin:$PATH"
treeknit test/NYdata/tree_ha.nwk test/NYdata/tree_na.nwk
OUT=$?
if [ "$OUT" != 0 ]; then
    rm -r treeknit_results
    echo CLI run fail
    exit 1
else
    rm -r treeknit_results
    echo CLI run OK
fi