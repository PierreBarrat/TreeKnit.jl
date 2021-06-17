# Overview

## Handling trees

Functions that direclty handle trees are found in the separate *TreeTools* package. 
  Here is a short list of useful ones: 
  - `read_tree(file)`: read tree from newick file. 
  - `parse_newick(string)`: parse newick `string` into a `TreeNode` object
  - `node2tree(n::TreeNode)`: create a `Tree` object from node `n`, using it as a root. 
  - `write_newick(file::String, t::Tree)`/`write_newick([file::String], n::TreeNode)`: write tree to `file` using newick format. Return a newick string if `file` is not provided. 

## Simple case 

Let's see how to infer MCCs for a very simple case: two trees with five leaves. 
```@example basic; continued = true 
using RecombTools
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

Note that the output of `computeMCCs` is a `Dict`, indexed by pairs of keys of the input dictionary `trees`. 

By convention, `mccs[i,i]` is always empty.

## More than two trees
If more than two trees are given as input, `computeMCCs` infers MCCs for all pairs of trees.  
```@example more_trees
using RecombTools # hide
t1 = node2tree(parse_newick("((A,B),((C,Y),(D,X)))"))
t2 = node2tree(parse_newick("((A,(B,X)),((C,Y),D))"))
t3 = node2tree(parse_newick("((A,(B,Y)),(C,(D,X)))"))
trees = Dict(1=>t1, 2=>t2, 3=>t3)
mccs = computeMCCs(trees)
mccs[1,2]
```
```@example more_trees
mccs[1,3]
```

## [Naive estimation](@id naive_mccs)
It is also possible to compute a "naive" estimation of MCCs using the `naive` keyword. 
  When `naive` is set to `true`, `computeMCCs` returns maximum clades that are exactly compatible between pairs of trees: 
```@example naive
using RecombTools # hide
t1 = node2tree(parse_newick("(((A1,A2),(B1,B2)),(C1,C2))"))
t2 = node2tree(parse_newick("(((A1,A2),(C1,C2)),(B1,B2))"))
trees = Dict(1=>t1, 2=>t2)
computeMCCs(trees; naive=true)
```