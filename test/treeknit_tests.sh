all_tests=0

julia --project=. deps/build.jl
retval="$?"
if [ "$retval" == 0 ]; then
    echo Installation OK
else
    echo Installation FAIL
    ((all_tests++))
fi

export PATH="$HOME/.julia/bin:$PATH"
treeknit test/NYdata/tree_ha.nwk test/NYdata/tree_na.nwk
retval="$?"
if [ "$retval" == 0 ]; then
    echo CLI run OK
else
    echo CLI run FAIL
    ((all_tests++))
fi
rm -r treeknit_results

julia --project=. test/runtests.jl 
retval="$?"
if [ "$retval" == 0 ]; then
    echo tests OK
else
    echo tests FAIL
    ((all_tests++))
fi

if [ "$all_tests" == 0 ];then
	echo "All tests passed"
	exit 0
else
	exit "$all_tests"
fi
