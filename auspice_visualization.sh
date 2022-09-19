#!/usr/bin/env bash -l

eval "$(conda shell.bash hook)"
if conda env list | grep ".*nextstrain.*" >/dev/null 2>&1; then
    conda activate nextstrain
else
    conda env create --name nextstrain
    conda activate nextstrain
    exit
fi;
conda install --yes nextstrain-cli

augur export v2 -t "$1.resolved.nwk" --node-data "branch_lengths_$1.json" --output "$1.json"
augur export v2 -t "$2.resolved.nwk" --node-data "branch_lengths_$2.json" --output "$2.json"
auspice view --datasetDir $PWD
