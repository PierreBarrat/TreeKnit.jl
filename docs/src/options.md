# [Options](@id options)

`computeMCCs` uses a topogy based heuristic optimization to find maximally compatible clades. 
  Options to this heuristic are provided through the `OptArgs` object, that is optionally passed as a second argument.
  Essential options are detailed here. 	

## [Parsimony parameter](@id gamma)
The heuristic method used by *TreeKnit* tries to prune consistent clades from a pair of trees in order to increase a compatibility score between other clades. 
  Pruning a clade is interpreted as fixing a reassortment right above it, while increasing the compatibility between remaining clades removes reassortments. 
  A purely parsimonious heuristic should thus give the same weight to fixing a reassortment through pruning a clade and fixing one incompatibility in the trees. 

Here, we assign the score $\gamma$ to each pruned clade, and count as $1$ each incompatibility fixed for the remaining clades. 
  For $\gamma=1$, we obtain a parsimonious method that attempts to minimize the overall number of reassortments. 
  For higher values, pruning a clade must "fix" at least $\gamma$ incompatibilities to be considered a good move, making the obtained MCCs less parsimonious. 
  For infinite $\gamma$, pruning clades is impossible, and we fall back on the [naive estimation of MCCs](@ref naive_mccs). 
  
The example below illustrates the difference between different $\gamma$ values: 
```@example gamma1
using TreeKnit # hide
t1 = node2tree(parse_newick("((((A,B),C),D),E)"))
t2 = node2tree(parse_newick("((((D,B),E),A),C)")) # Same topology, but shuffled leaves
nothing # hide
```
Here, pruning the two leaves `(A,C)` or `(D,E)` results in compatible trees (resp. `((B,D),E)` and `((A,B),C)`). 
  These moves each have a cost 2$\gamma$ (removing 2 clades), but bring us from trees with 5 incompatibilities to 0. 
  They will only be accepted if $\gamma \leq 2.5$. 

```@repl gamma1
  computeMCCs(t1, t2, OptArgs(γ=2))
  computeMCCs(t1, t2, OptArgs(γ=3))
```


## Resolving trees with polytomies 
See [Resolving](@ref resolving)

## [Degeneracy: sorting with likelihood](@id likelihood)
When several MCC decompositions are possible, degeneracy is removed by using the `likelihood_sort` option (activated by default). 
In the example below, there are three equivalent decompositions if only topology is considered: 
```@example degeneracy
using TreeKnit # hide
t1 = node2tree(parse_newick("((A:2,B:2):2,C:4)"))
t2 = node2tree(parse_newick("(A:2,(B:1,C:1):1)"))
oa = OptArgs(likelihood_sort = false)
unique([computeMCCs(t1, t2, oa) for rep in 1:50]) # Repeating computation many times 
```

When taking branch lengths into account, this degeneracy vanishes: 
```@example degeneracy
oa = OptArgs(likelihood_sort = true)
unique([computeMCCs(t1, t2, oa) for rep in 1:50])
```