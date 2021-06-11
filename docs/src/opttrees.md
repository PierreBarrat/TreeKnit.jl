# The `opttrees` function

The core of the heuristic *RecombTools* is based on happens in the `opttrees` function, found in the `SplitGraph` submodule. 
  Given two trees, `opttrees` attempts to reconcile them by pruning certain clades. 
  A quick description of how the function operates is given here, with the two simple trees below as an example case: 
```@example opttrees
using RecombTools# hide
nwk1 = "(((A1:1,A2:2):2,(B1:2,(B2:1,B3:1):1):2):2,(C1:1,C2:2):4)";
nwk2 = "((A1:1,A2:2):2,((B1:2,(B2:1,B3:1):1):1,(C1:1,C2:2):1):1)";
t1 = node2tree(parse_newick(nwk1))
t2 = node2tree(parse_newick(nwk2))
nothing # hide
```
Trees are not displayed here for space reasons, but you're encouraged to draw them if you want to follow along! 

## Coarse-graining of naive MCCs

As a first step, naive MCCs are computed for input trees using the `naive_mccs` function. 
  Here, we find three clades that are already compatible between the two trees: 
```@example opttrees
		treelist = Any[t1, t2]
		mcc = naive_mccs(treelist)
```  

Trees are then "reduced" to those MCCs: new trees are built where each leaf corresponds to one of the naive MCCs. 
  The reduced trees have incompatibilities at the leaf level: it is no longer possible to group some of their leaves together in a consistent clade.  
```@example opttrees
mcc_names = RecombTools.name_mcc_clades!(treelist, mcc)
for (i,t) in enumerate(treelist)
	treelist[i] = RecombTools.reduce_to_mcc(t, mcc)
end
nothing # hide
```

The trees in `treelist` are now a reduced form of `t1` and `t2`, and the names of the new leaves correspond to clades in the original tree. 
  The mapping between leaf name and original clade is stored in `mcc_names`
```@repl opttrees
treelist[1]
treelist[2]
mcc_names
```

## `SplitGraph` object and optimization

Once the trees reduced to their naive MCCs, we construct a `SplitGraph` object from them. 
- each leaf indexed by an integer
- each internal node stores the ensemble of leaves that are underneath it
- internal nodes have a unique integer "color" that corresponds to the index of the tree they belong to, *e.g.* 1 for `t1` and 2 for `t2`. 
