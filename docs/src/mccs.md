# Maximally Compatible Clades

*TreeKnit* reconstructs the ARG from trees by first inferring Maximally Compatible Clades (MCC). 

## Handling trees

Functions that directly handle trees are found in the separate *TreeTools* package. 
  Here is a short list of useful ones: 
  - `read_tree(file)`: read tree from newick file, return `Tree` object. 
  - `parse_newick(string)`: parse newick `string` into a `TreeNode` object.
  - `node2tree(n::TreeNode)`: create a `Tree` object from node `n`, using it as a root. 
  - `write_newick(file::String, t::Tree)`/`write_newick([file::String], n::TreeNode)`: write tree to `file` using newick format. Return a newick string if `file` is not provided. 

## Simple case 

Let's see how to infer Maximally Compatible Clades (MCC) for a very simple case: two trees with five leaves. 
```@example basic; continued = true 
using TreeKnit
t1 = node2tree(parse_newick("((A,B),(C,(D,X)))"))
t2 = node2tree(parse_newick("((A,(B,X)),(C,D))"))
```

The `computeMCCs` function takes two trees as input. 
```@example basic
mccs = computeMCCs(t1, t2)
```
Individual MCCs are simply arrays containing labels of leaves of the trees.  

## Interpretation of results

The genealogy of two RNA segments subject to reassortment is described by an Ancestral Reassortment Graph (ARG). 
An ARG is a directed graph that represents the lineage of a given pair of segments by coalescence of nodes, as in a genealogical tree, but also shows reassortment events and the exchange of segments by nodes that have two ancestors. 
Since reassortments only occur between segments, the genealogy of given segment is described by a tree. 
As a result, the ARG must embed both segment-trees, and every branch in the ARG has to belong to either one of the trees, or to both. 

*TreeKnit* infers the ARG by finding the branches that are common to both trees. 
Given two trees with potentially different topologies, it tries to "glue" them together in a reasonable way, where the interpretation of reasonable can vary between *parsimonious* and *conservative* (see the [parsimony parameter](@ref gamma) $\gamma$). 

The MCCs returned by `computeMCCs` represent regions of the ARG (and of the segment trees) where branches are common to both trees. 
In other words, these are the regions where the two segment trees must be "glued together". 
Given those regions and the knowledge of the trees, it is possible to unambiguously reconstruct the genealogy. 

!!! info "Number of reassortments in the genealogy"
    When going up the ARG (backwards in time), a reassortment consists of passing from a region where branches are common to the two trees to a region where they are not. It is a *split* of branches. 
    As a consequence, the root of each MCC must be a reassortment, *with the exception* of an MCC containing the root of both trees. 
    The number of reassortments events in the inferred ARG can thus simply be obtained by counting the number of MCCs, potentially removing the one that contains the roots of both trees if it exists. 




## More than two trees
If more than two trees are given as input, `computeMCCs` infers MCCs for all pairs of trees.  
Use a `Dict` to pass the trees as input. 
```@example more_trees
using TreeKnit # hide
t1 = node2tree(parse_newick("((A,B),((C,Y),(D,X)))"))
t2 = node2tree(parse_newick("((A,(B,X)),((C,Y),D))"))
t3 = node2tree(parse_newick("((A,(B,Y)),(C,(D,X)))"))
trees = Dict(1=>t1, 2=>t2, 3=>t3)
mccs = computeMCCs(trees)
mccs[1,2] # MCCs for t1 and t2
```
```@example more_trees
mccs[1,3] # MCCs for t1 and t3
```

The output `mccs` is also a `Dict`, indexed by pairs of keys of the input dictionary `trees`. 
Indexing is symmetric: `mccs[i,j] == mccs[j,i]`. 
By convention, `mccs[i,i]` exists and is empty. 

## [Naive estimation](@id naive_mccs)
It is also possible to compute a "naive" estimation of MCCs using the `naive` keyword. 
  When `naive` is set to `true`, `computeMCCs` returns maximum clades that are exactly compatible between pairs of trees: 
```@example naive
using TreeKnit # hide
t1 = node2tree(parse_newick("(((A1,A2),(B1,B2)),(C1,C2))"))
t2 = node2tree(parse_newick("(((A1,A2),(C1,C2)),(B1,B2))"))
trees = Dict(1=>t1, 2=>t2)
computeMCCs(trees; naive=true)
```