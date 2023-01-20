# [Resolving](@id resolving)

The lack of resolution in the inference of phylogenetic trees results in *polytomies*: internal nodes with more than two offsprings. 
  Polytomies can cause the topology of two trees to differ, which cause problems when inferring reassortments using topological information. 
  Suppose for instance that we have a first tree (*e.g.* for a given segment of flu): 
```@example 1
using TreeKnit # hide
t1 = parse_newick_string("(A,(B,C));")
```
where we can identify the clade `(B,C)` because of a mutation in this segment present in `B` and `C` but not in `A`. 

If this mutation does not exist in a second segment, then in the absence of reassortment its tree will look something like this: 
```@example 1
t2 = parse_newick_string("(A,B,C);")
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
t1 = parse_newick_string("(A,(B,C));")
t2 = parse_newick_string("(A,B,C);")
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

If the `resolve` option is passed to `run_treeknit!` (through `OptArgs`), `resolve!` will be called on each pair of trees before each iteration of the MCC inference procedure, MCC inference will allow for resolution and after each pair-wise iteration of `TreeKnit` the trees shall be resolved using the inferred MCCs.



## Resolving during MCC inference
The above resolving method works fine for simple "obvious" cases, where a split in one tree directly resolves a polytomy in another. 
However, consider the following case: 
```@setup 3
using TreeKnit
```
```@example 3; continued = true
t1 = parse_newick_string("((A,B),(C,(D,(E,X))));")
t2 = parse_newick_string("((A,(B,X)),(C,D,E));")
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
run_treeknit!(t1, t2, OptArgs(;resolve=false))
```

In order to achieve progress in this kind of situation, we have to perform two operations at the same time: 
  - realize that `X` is the only reassorted strain, and can be ignored when resolving.
  - resolve `t2` with the `(D,E,X)` split, ignoring `X`. 

This is done automatically during MCC inference if the `resolve` option of `OptArgs` is given (default):  
```@example 3
MCCs = run_treeknit!(t1, t2, OptArgs(;resolve=true))
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

## [Strict vs liberal resolve](@id resolve_strict_vs_liberal)

When resolving using MCCs, TreeKnit uses one of two options: strict (default) or liberal resolution (`--liberal-resolve` flag). 
Given two trees `t1` and `t2` and their MCCs, it is not always possible to unambiguously resolve `t2` using the splits of `t1` *even* in an MCC. 
To explain this, we consider the two following trees: 
```@repl strict_lib
using TreeKnit # hide
using TreeTools # hide
t1 = parse_newick_string("((A,(B,C)),D);")
t2 = parse_newick_string("(A,B,C,D);")
```

Assume that the MCCs are `[A,B,C]` and `[D]`. In this minimal example this is a bit contrived, but this type of situation can take place in larger  trees. 
The essential ingredients here are 
- the polytomy `(A,B,C,D)` in `t2`. 
- the fact that the `(A,(B,C))` clade is nicely resolved in `t1`
- `(A,B,C)`, `D`, and the other nodes being in different MCCs, meaning there is a reassortment above leaf `D` and above the MRCA of `(A,B,C)`

At first sight, we would like to introduce the splits `(A,B,C)` and `(B,C)` into `t2`: these splits exist in `t1` and are in the shared region `(A,B,C)`. 
This would result in the following for `t2`: 

```@example strict_lib
t2_resolved_1 = parse_newick_string("((A,(B,C)),D);")
```

However, since there is a reassortment above `D`, we cannot exclude `D` being nested in the `(A,B,C)` clade 
This gives us other possibilities for `t2` (non-exhaustive):

```@repl strict_lib
t2_resolved_2 = parse_newick_string("(A,(B,(C,D)));")
t2_resolved_3 = parse_newick_string("(A,((B,C),D));")
```

All the above possibilities to resolve `t2` are compatible with the found MCCs and the splits in `t1`. 

The `strict` and `liberal` options for resolution make different choices in this situation. 
- Strict resolve will consider this situation as ambiguous, and not attempt any resolution. Hence the final `t2` will have a Newick string `"(A,B,C,D);"`. 
  More generally, with the `strict` option, TreeKnit ill only resolve a polytomy if the relation of all branches within the polytomy can be unambiguously determined using the other tree, potentially not fully resolving shared regions of the trees. 

```@repl strict_lib
MCCs = [["D"], ["A", "B", "C"]]
new_splits_strict = resolve!(t1, t2, MCCs; strict=true);
isempty(new_splits_strict)
t2
```

- Liberal resolve will resolve trees as much as possible, fully resolving shared regions of the two trees, but arbitrarily choosing the location of these aforementioned branches, leading to potentially wrong splits. 
  As in the example below, liberal resolve will always try to "pull" MCCs out of polytomies 
```@example strict_lib
t1 = parse_newick_string("((A,(B,C)),D);") # hide
t2 = parse_newick_string("(A,B,C,D);") # hide
new_splits_liberal = resolve!(t1, t2, MCCs; strict=false);
new_splits_liberal[2]
t2
```

Here are a few extra notes on strict vs liberal options: 
- liberal resolve will in general result in significantly more wrong splits placed in the trees
- to create an ARG, it is necessary that MCCs are resolved exactly the same in both trees: the trees must be *knitted* together in those regions. 
  For this reason, it is also necessary to use liberal resolve in this case. 