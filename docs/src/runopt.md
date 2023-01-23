# The `runopt` function

```@setup runopt
using TreeKnit
using TreeKnit.SplitGraph
using TreeTools
```

For larger trees, `opttrees` will not find all reassortments at once. 
  The reason for this is that it is often necessary to "clean" reassorted leaves in smaller clades in order to see other reassortments deeper in the tree. 
  As a consequence, trees obtained after one round of `opttrees` will still have incompatibilities. 
Below is a relatively simple example of trees for which `opttrees` does not find all reassorted leaves in one go: 

```@example runopt
nwk1 = "(Z,(G,(((A,X),(B,C)),((D,Y),(E,F)))));"
nwk2 = "(G,((A,(B,(C,X))),((D,(E,(F,Y))),Z)));"
t1 = parse_newick_string(nwk1; label="t1")
t2 = parse_newick_string(nwk2; label="t2")
nothing # hide
```

These trees are constructed in the following way (you're encouraged to draw them!): 
- Clade `((A,X),(B,C))` in the first tree "corresponds" to clade `(A,(B,(C,X)))` in the second tree, with `X` as the reassorted leaf. 
  Same goes for clades `((D,Y),(E,F))` and `(D,(E,(F,Y)))`. 
- We respectively name these not yet compatible clades `ABC` and `DEF`. At a deeper level, the trees are now of the form `(Z,(G,(ABC,DEF)))` for the first and `(G,(ABC,(DEF,Z)))` for the second, with `Z` as the obvious reassorted leaf. 

The important property here is that for a high enough value of $\gamma$, it is only possible for `opttrees` to see that `Z` is reassorted if clades `ABC` and `DEF` are "coarse-grained" to leaves. 
  In return, this coarse-graining is only possible after `X` and `Y` have been identified as reassorted leaves, which will happen after a first iteration of `opttrees`.  
  The task of iterating `opttrees` is performed by the `runopt` function. 
  Here, we walk through typical steps it takes.


Let's do a first pass with `opttrees`, with $\gamma=3$.  
  We first build the `SplitGraph` object, as detailed in the [`opttrees`](@ref opttrees) page. 

```@example runopt
treelist = Any[t1, t2]
mccs_naive = naive_mccs(treelist) # these are just the leaves in this example
mcc_names = TreeKnit.name_mcc_clades!(treelist, mccs_naive)
for (i,t) in enumerate(treelist)
  treelist[i] = TreeKnit.reduce_to_mcc(t, mccs_naive)
end
g = SplitGraph.trees2graph(treelist);
nothing # hide
```

We then run the simulated annealing optimization to find optimal leaves to remove. 

```@example runopt
opt_confs = SplitGraph.sa_opt(g; Trange = reverse(1e-3:1e-2:1), M = 10, γ = 3)[1]
mccs_found = [mcc_names[x] for x in g.labels[.!opt_confs[1]]]
```

!!! info $\gamma = 2$
    If you run this example using $\gamma \leq 2$, `Z` will immediatly be found as a reassorted strain. 
    It is indeed not easy to find example that combine trees with a small number of leaves (9 here), obvious reassortments, and that are not solved in one go by `opttrees` with a low value of $\gamma$. 
    However, when dealing with trees with hundreds of leaves, finding all MCCs in one go is the exception rather than the rule. 


As expected, `X` `Y` are found as reassortants. 
  However, the two trees will still have incompatibilities when removing those two leaves. 
  To make this explicit, we remove the leaves `X` and `Y` and compute naive mccs again. 

```@example runopt
TreeKnit.pruneconf!(mccs_found, treelist...) # prune clades in a list of trees. Wrapper around TreeTools.prunesutree!
mccs_naive = naive_mccs(treelist...)
```

Now that `X` and `Y` are removed, we see that clades `ABC` and `DEF` are common to both trees. 
  If we reduce the pruned trees to their new naive MCCs again, we now see that `Z` is an obvious choice for a reassorted strain: 

```@repl runopt
mcc_names = TreeKnit.name_mcc_clades!(treelist, mccs_naive)
for (i,t) in enumerate(treelist)
  treelist[i] = TreeKnit.reduce_to_mcc(t, mccs_naive)
end
treelist[1]
treelist[2]
``` 

To finish the inference of MCCs, we would now have to re-run the optimization process. 
  It is now clear that the `opttrees` has to be iterated. 
  This is performed automatically by the `runopt` function. 
  This process stops when one of the following end conditions is found: 

1. If no new MCCs are found in a given iteration of `opttrees`. 
   This occurs when the optimal configuration resulting from the simulated annealing has all leaves present. 
2. If not and new MCCs were found, prune them from the trees. 
   If the resulting trees do not have any incompatibility.
3. If not, if the maximum number of iterations has been reached. 
   This can be set through `OptArgs`, with a default of 15. 

 