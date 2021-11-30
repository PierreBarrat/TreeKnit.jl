# [The `opttrees` function](@id opttrees)

The core of the heuristic *TreeKnit* is based on happens in the `opttrees` function, found in the `SplitGraph` submodule. 
  Given two trees, `opttrees` attempts to reconcile them by pruning certain clades. 
  A quick description of different steps in this function is given here, with the two simple trees below as an example case: 
```@example opttrees
using TreeKnit# hide
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
mcc_names = TreeKnit.name_mcc_clades!(treelist, mcc)
for (i,t) in enumerate(treelist)
	treelist[i] = TreeKnit.reduce_to_mcc(t, mcc)
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

## The `SplitGraph` object 

Once the trees reduced to their naive MCCs, we construct a `SplitGraph` object from them. 

!!! info "*SplitGraph* submodule"
    The `SplitGraph` type and some of the functions used below are in the *SplitGraph* submodule of *TreeKnit*. Access them by calling `using TreeKnit.SplitGraph` and preceding the calls by `SplitGraph.`

The `SplitGraph` is a directed graph that is based on both trees, and has two kind of nodes: 
- leaf nodes correspond to leaves of the trees, and are identifier by integers. 
  They have as many ancestors as there are trees in the `SplitGraph`. 
- internal nodes, called `SplitNode`s, correspond to internal nodes in one of the two trees. A `color::Int` attribute identifies the tree to which they belong (*e.g.* 1 for the first tree, 2 for the second, etc...).
  They have only one ancestor, of the same color. 
  Importantly, they are identified by the ensemble of leaf nodes that are "below" them, that is the subset of all their direct and indirect offsprings that are leaves. 
  As such, they uniquely correspond to a split in one of the two trees. 
  This information is stored as an array of integer in their `conf` field. 

Let us know build the `SplitGraph` object: 
```@example opttrees
using TreeKnit.SplitGraph
g = SplitGraph.trees2graph(treelist); 
g.labels_to_int
```

!!! warning
    It is recommemded that you add `;` to the end of lines when working with `SplitGraph`, `SplitNode` or `LeafNode` in the REPL. 
    If you forget, you will quickly see why :-) 

We now see that the three leaves from our coarse-grained trees have been attributed an integer index in the `SplitGraph`. 
  Let us take a look at the internal nodes above the leaf `MCC_3`: 
```@repl opttrees
a1 = g.leaves[g.labels_to_int["MCC_3"]].anc[1]; # Ancestor for the first tree
a2 = g.leaves[g.labels_to_int["MCC_3"]].anc[2]; # Ancestor for the second tree
[a1.color, a2.color] # a1 and a2 resp. belong to trees 1 and 2
a1.conf # list of leaves below `a1`. Among those is the index for "MCC_3".
[g.labels[i] for i in a1.conf] # Same as above, with labels
[g.labels[i] for i in a2.conf] # and the same for a2 
```
We now immediatly see that the internal nodes above `MCC_3` in the two trees define different splits: `(MCC_1, MCC_3)` in the first tree is different from `(MCC_2, MCC_3)` in the second tree. 
  This is the idea underlying the inference of MCCs. 

## Counting incompatibilities

In the example above, the ancestors of leaf `MCC_3` in the two trees define different splits: this is called an incompatibility. 
  Examination of the trees reveals that there are also similar incompatibilitie for the two other leaves `MCC_1` and `MCC_2`. 
  This can be computed using the `count_mismatches` function: 

```@example opttrees
SplitGraph.count_mismatches(g)
```

We indeed find 3 mismatches, one for each leaf. 

However, it is possible to explain the two example trees using less than three reassortments. 
  What would happen for example if we removed `MCC_1` from both trees? 
  The first non-trivial split above leaf `MCC_3` in both trees would then be `(MCC_2, MCC_3)`, and the same goes for leaf `MCC_2`. 
  The number of incompatibilities would then go down to 0. 
"Removing" leaves from the trees, or the graph, is done by defining a *configuration*: an array of booleans that stores the presence or absence of each leaf. 
  To remove `MCC_1`, we simply design a configuration that has `0` at the index corresponding to `MCC_1`: 

```@example opttrees
  conf = ones(Bool, length(g.leaves))
  conf[g.labels_to_int["MCC_1"]] = false # Remove `MCC_3` from the configuration
  conf
```

To compute the number of incompatibilities given a configuration, we use the `compute_energy` function. The result is interpreted as the "energy" of this configuration given the graph `g`: 

```@example opttrees
SplitGraph.compute_energy(conf, g)
```

!!! info 
    The function `count_mismatches(g)` shown above is a simple shortcut for 
    ```julia
    conf = ones(Bool, length(g.leaves))
    SplitGraph.compute_energy(conf, g)
    ```
    In other words, it computes the energy for the configuration where all leaves are present. 

By removing a leaf, *i.e.* by "enforcing" a reassortment right above it, we've reduced the number of incompatibilities for the remaining ones to 0. 
  Since removing a leaf corresponds to "enforcing" a reassortment, we have to assign a cost to it, that we call $\gamma$. 
  This defines a score for each configuration, defined as the difference between the energy of the configuration and $\gamma$ times the number of leaves that were removed.
Depending on the value of $\gamma$, the difference in overall score associated to removing a leaf or keeping it will change from negative to positive. 
  Scores are computed with the `compute_F` function that takes $\gamma$ as its last argument.
  Here are the differences in scores before and after removing `MCC_3`, for different values of $\gamma$:   
```@repl opttrees
conf0 = ones(Bool, length(g.leaves)) # Configuration with all leaves
SplitGraph.compute_F(conf, g, 1) - SplitGraph.compute_F(conf0, g, 1)
SplitGraph.compute_F(conf, g, 2) - SplitGraph.compute_F(conf0, g, 2)
SplitGraph.compute_F(conf, g, 3) - SplitGraph.compute_F(conf0, g, 3)
SplitGraph.compute_F(conf, g, 4) - SplitGraph.compute_F(conf0, g, 4)
```

For this simple example, $\gamma = 3$ is the "critical" value above which the fact of removing `MCC_3` or any other leaf is not considered a good move. 
  The inference of MCCs for $\gamma \leq 3$ and $\gamma > 3$ will thus give different results. 
  In the first case, two MCCs will be found, corresponding to one reassortment event (above `MCC_3` for instance)
. 
  In the second, three MCCs and three reassortments will be found. 

```@repl opttrees
trees = Dict(1=>t1, 2=>t2);
computeMCCs(trees, OptArgs(γ=3.1))[1,2]
computeMCCs(trees, OptArgs(γ=2.9))[1,2]
```

## Simulated annealing 

The `opttrees` function attempts to find the configuration, *i.e.* a set of leaves to remove, that minimizes the compatibility score presented above. 
  Since this is a discrete optimization problem with no clear mathematical formalization, we choose to use the [simulated annealing](https://en.wikipedia.org/wiki/Simulated_annealing) technique. 
Let us find optimum configurations for our simple trees: 

```@example opttrees
Trange = reverse(1e-3:1e-2:1) # Cooling schedule
M = 10 # Number of iterations per temperature value
opt_confs = SplitGraph.sa_opt(g; Trange, M, γ = 2)[1]
```

We find three optimal configurations, each corresponding to removing one leaf. 
  This indeed corresponds to the three possible single-reassortment explanations that we could give to reconcile the two trees. 
  Without branch length information, it is impossible to choose between one of these three optimas: the problem is degenerate. 
  A likelihood based way to break this degeneracy using branch length is described [here](@ref likelihood).

For now, let's imagine that we have chosen the first optimum configuration as our best solution. 
  Let's now map it back on the initial trees: 

```@repl opttrees
removed_leaves = g.labels[.!opt_confs[1]] # Expressed with coarse grained leaves
# `mcc_names` was defined above
removed_clades = [mcc_names[x] for x in removed_leaves]
```

What this means is that we have just inferred all elements in `removed_clades` (just one in our case) to be MCCs. 
  Of course, in this simple example, it is immediate to see that the other MCC simply consists of all the remaining leaves in the original trees. 
  This can also be deduced from the fact that the energy of all the optimal configurations is 0. 
However, in the general case, some incompatibilities will remain even after simulated annealing. 
  For this reason, the `opttrees` function only outputs MCCs that have been identified by having removed them from the tree, *i.e.* by having enforced a reassortment above their root node. 
  If the trees that remain after having pruned these MCCs still have incompatibilities, the process described here needs to be *iterated*. 
  This is performed by the `runopt` function. 






