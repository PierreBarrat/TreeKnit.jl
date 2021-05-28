# Usage

## Simple case 

Let's see how to infer MCCs for a very simple case: two trees with five leaves. 
```@example basic; continued = true 
using RecombTools, TreeTools
t1 = node2tree(parse_newick("((A,B),(C,(D,X)))"))
t2 = node2tree(parse_newick("((A,(B,X)),(C,D))"))
```

The `computeMCCs` function takes a dictionary of trees as input. 
It would normally be indexed by *e.g.* flu segments, but here we will simply index it with integers.

```@example basic
trees = Dict(1=>t1, 2=>t2)
mccs = computeMCCs(trees)
mccs[1,2]
```

Individual MCCs are simply arrays containing labels of leaves of the trees.  

Note that the output of `computeMCCs` is a dictionary, indexed by pairs of keys of the input `trees`. 
This is mostly practical when computing MCCs for more than two trees: 

```@example basic
t3 = node2tree(parse_newick("((A,B),(C,(D,X)))")) # t3 is equal to t1
trees[3] = t3
mccs = computeMCCs(trees)
mccs[1,3]
```

By convention, `mccs[1,1]` and `mccs[2,2]` are empty. 


## Degeneracy
```@example degeneracy
using RecombTools, TreeTools # hide
t1 = node2tree(parse_newick("((A:2,B:2):2,C:4)"))
t2 = node2tree(parse_newick("(A:2,(B:1,C:1):1)"))
trees = Dict(1=>t1, 2=>t2)
oa = OptArgs(likelihood_sort = false)
unique([computeMCCs(trees, oa)[1,2] for rep in 1:10])
```

```@example degeneracy
oa = OptArgs(likelihood_sort = true)
unique([computeMCCs(trees, oa)[1,2] for rep in 1:10])
```