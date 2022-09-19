julia --project=. deps/build.jl
if [ $? -eq 0 ]; then
    echo Installation OK
else
    echo Installation FAIL
fi

export PATH="$HOME/.julia/bin:$PATH"
treeknit test/NYdata/tree_ha.nwk test/NYdata/tree_na.nwk
if [ $? -eq 0 ]; then
    echo CLI run OK
else
    echo CLI run FAIL
fi

julia --project=. test/runtests.jl 
if [ $? -eq 0 ]; then
    echo tests OK
else
    echo tests FAIL
fi
