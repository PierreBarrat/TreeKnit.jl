# Resolving (@id resolving)

The lack of resolution in the inference of phylogenetic trees results in *polytomies*: internal nodes with more than two offsprings. 
  Polytomies can cause the topology of two trees to differ, which cause problems when inferring reassortments using topological information. 
  Suppose for instance that we have a first tree (*e.g.* for a given segment of flu): 
```@example 1
using TreeKnit # hide
t1 = node2tree(parse_newick("(A,(B,C))"))
```
where we can identify the clade `(B,C)` because of a mutation in this segment present in `B` and `C` but not in `A`. 

If this mutation does not exist in a second segment, then in the absence of reassortment its tree will look something like this: 
```@example 1
t2 = node2tree(parse_newick("(A,B,C)"))
```
As a result, `t1` and `t2` differ for a reason unrelated with reassortment. 

To overcome this issue, one has to **resolve** trees as much as possible, typically using the information of one to remove polytomies in the other. 
  This can only be done if no reassortments are present. 

## Resolving pairs of trees

In the example above, it is natural to think that the difference between the two trees is due to a lack of resolution and not to reassortment. 
  This is because the split `(B,C)` in the first tree is **compatible** with the second tree: it is possible to add this split to `t2`. 
  What this means is that the difference in topology can be explained by something else than reassortment. 

In this case, we can simply resolve `t2` by adding each split in `t1` with which it is compatible. 
  If `t1` has polytomies, the same could be done to resolve `t1` using `t2`. 
  This operation is performed by the `resolve!` function: 
```@setup 2
using TreeKnit
t1 = node2tree(parse_newick("(A,(B,C))"))
t2 = node2tree(parse_newick("(A,B,C)"))
```
```@example 2
new_splits = resolve!(t1, t2);
t2
```

`resolve!` returns an array of `SplitList` objects (of the *TreeTools* package) containing the new splits introduced in each of the two trees: 
```@repl 2
isempty(new_splits[1])
new_splits[2]
```

If the `resolve` option is passed to `computeMCCs` (through `OptArgs`), `resolve!` will be called on each pair of trees before each iteration of the MCC inference procedure. 

## Resolving during MCC inference
The above resolving method works fine for simple "obvious" cases, where a split in one tree directly resolves a polytomy in another. 
However, consider the following case: 
```@setup 3
using TreeKnit
```
```@example 3; continued = true
t1 = node2tree(parse_newick("((A,B),(C,(D,(E,X))))"))
t2 = node2tree(parse_newick("((A,(B,X)),(C,D,E))"))
```
There are now two sources of topological differences between `t1` and `t2`: 
- The reassorted strain `X`. 
- The lack of resolution resulted in a polytomy `(C,D,E)` in `t2`, which is resolved in `t1` in the form of `(C,(D,(E,X))` 

The `resolve!` function is helpless in such cases: 
```@example 3
new_splits = resolve!(t1, t2)
isempty(new_splits[2])
```
Indeed, there is no split in `t1` that can directly help us resolve `t2`. 
  The closest such split is `(D,E,X)`, but it is incompatible with `t2` because of `X`. 
  If we knew beforehand that `X` is reassorted, we could simply ignore it while resolving `t2`. 
  The `(D,E,X)` split in `t1` would become `(D,E)`, which is compatible with `t2`, and the `resolve!` function would handle this.  

However, the topology-based heuristic used by *TreeKnit* is not able to detect that `X` is the only reassorted leaf *if the trees are not resolved*!
  Indeed, if we "remove" `X` from the trees, some incompatibilities will remain. 
  For instance, the split above `E` will be `(D,E)` in the first tree and `(C,D,E)` in the second. 
  Without resolving, the heuristic will predict a reassortment above almost every leaf: 
```@example 3
computeMCCs(t1, t2, OptArgs(;resolve=false))
```

In order to achieve progress in this kind of situation, we have to perform two operations at the same time: 
  - realize that `X` is the only reassorted strain, and can be ignored when resolving.
  - resolve `t2` with the `(D,E,X)` split, ignoring `X`. 

This is done automatically during MCC inference if the `resolve` option of `OptArgs` is given (default):  
```@example 3
MCCs = computeMCCs(t1, t2, OptArgs(;resolve=true))
```

## Resolving with inferred MCCs

Once the MCCs are inferred, it is possible to use them to resolve trees: in the regions of shared branches of the ARG, the two trees `t1` and `t2` must have the same splits. 
  The `resolve!` function also has a method for this. 
  Using the example above, we have
```@repl 3
resolved_splits = resolve!(t1, t2, MCCs)
t2
```

The split `(D,E)` is now present in `t2`. 
Note that it was not present in `t1`: only the splits `(D,E,X)` and `(E,X)` existed there. 
However, since `resolve!` now knows `(A,B,C,D,E)` is an MCC, the `resolve!` function can "ignore" leaf `X` when resolving.   