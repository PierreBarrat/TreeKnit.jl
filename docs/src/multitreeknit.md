# [MultiTreeKnit] (@id multitreeknit)

`TreeKnit` can be used to infer reassortment events between multiple tree pairs. When reassortment events between three or more trees should be inferred, `TreeKnit` uses a recursive inference strategy to:
- resolve each tree using information from all other trees in a consistent manner,
- identify shared regions of each tree pair (MCCs),
- find reassortment events between tree pairs.
The following example illustrates the benefits of this inference strategy over running `TreeKnit` individually on all tree pairs. Note however, that MultiTreeKnit may still return MCCs that are inconsistent with each other, which prevents the construction of an ARG. Therefore, we do not reconstruct an ARG for more than two trees. 

Furthermore, MultiTreeKnit has been optimized for using multiple reassorting segments to infer parameters of pathogen evolution. We optimized for sensitivity in MCC inference and polytomy resolution, minimizing the false positive rate for shared branches and new splits, but increasing the false negative rate, i.e. typically inferring too many reassortment events. 

We describe two versions of MultiTreeKnit, differing by the number of inference rounds. The first, default version uses only one round of inference and is optimal for tree reconstruction. The second, uses two rounds of inference, resolves more branches and produces slightly less accurate trees, but optimal MCCs.

Both approaches follow the same structure:

1. All trees are compatibly pre-resolved with each other by default (this can be deactivated with the `--no-pre-resolve` flag). This means that we compute each tree's set of splits and try to add each split from each tree to all other trees, if and only if that split is compatible with all other trees. This is similar to resolving trees prior to inference in the standard, two-tree `TreeKnit`.
2. Run `1` or `2` rounds of standard `TreeKnit` on all tree pairs. In the final round do not resolve the trees any further. Otherwise, if `2` rounds of `TreeKnit` are run, resolve tree pairs prior and during pair-wise inference in the first round, using `strict` resolve to only add unambiguous splits. If trees are resolved during inference tree pairs must be resolved in a specific order to avoid inconsistent resolution. To prevent MCCs calculated for tree pairs at the start of the round from no longer being accurate when trees are later resolved, the MCCs calculated in the initial round are discarded and a final round of pairwise `TreeKnit` is run on all resolved tree pairs.
3. Ladderize the first tree and sort all polytomies according to to this tree to allow for visualization as a tanglegram. 

![plot](./Pictures/MultiTK_example.png)
## Consistent Resolution

When `TreeKnit` is run individually on all tree pairs trees could be resolved inconsistently. By default `TreeKnit` resolves trees when searching for reassortment events, this is especially important for influenza where resolution is often low and not resolving trees can lead to much higher rates of reassortment being inferred. 

For the example above this would mean that when the MCCs of `tree a` and `tree c` are computed no reassortment events would be found as `tree c` would be resolved according to `tree a` (i.e. a branch would be introduced above leaves `A` and `B`). However, when the MCCs of `tree b` and `tree c` are computed `tree c` would be resolved according to `tree b` and a branch would be introduced above leaves `B` and `C`. These two different resolutions of `tree c` are incompatible with each other and do not allow us to use this reassortment event information together.

In our default version of `MultiTreeKnit` we would first try to pre-resolve all trees using each other. In this case, as the $(A, B)$ and the $(B, C)$ split are not compatible, we do not know which to introduce into `tree c` and thus we would not introduce any new splits. We would then compute all MCC-pairs without resolution, which would mean inferring reassortment events between all trees. 

The other approach would be to choose to resolve `tree c` with `tree a` OR with `tree b`, inferring less reassortment events, which is more likely to be true. However, we have no information which tree we should used to resolve `tree c`. By default, if we tell `MultiTreeKnit` to perform `2+` rounds of inference it will further resolve trees using tree order to choose how to resolve trees. The order that pairs are resolved in is shown in the picture below 

![plot](./Pictures/resolution_order.png)

Here $T_{1,2}$ corresponds to resolving `tree 1` and `tree 2` using each other and then generating a list of MCCs using standard `TreeKnit`. These resolved trees are then used to calculate the MCCs of the next neighboring tree pairs (which are connected by arrows). This means that after calculating the MCCs for `tree 1` and `tree 2` and resolving their polytomies using each other, instead of using the original `tree 1` the `resolved tree 1` is then further used when calculating the MCCs between `tree 1` and `tree 3`. Leading to a consistent tree resolution.

Furthermore, we use [strict resolution](@ref resolve_strict_vs_liberal) when resolving trees with MultiTreeKnit. Strict resolve will only resolve polytomies if the location of each branch in the new split can be fully determined by the other tree. Using strict resolve prevents the introduction of incorrect splits into the trees, this is especially important when resolved trees are used downstream for inference as these splits could prevent the simulated annealing from converging. 

At the end of the sequential inference on all tree pairs, each tree will be resolved as much as possible using each other tree. However, the output MCCs might not be consistent with each other and might not fulfill the necessary transitivity requirements to create an ARG. Furthermore, it can occur that the MCCs inferred for the tree pairs that were calculated at the start of the round are no longer topologically compatible now that those trees have been further resolved. Meaning that the inferred MCCs might not actually be MCCs. We prevent this from happening, and make sure that all MCCs that are inferred are actual MCCs by running a final round where we re-infer all tree pair MCCs without resolving trees. However, the MCCs might still not fulfill all necessary consistency conditions to produce an ARG.

## Consistent MCCs

When `TreeKnit` is run individually on all tree pairs the final MCCs might be inconsistent with each other. 
As can be seen in the example when `TreeKnit` is run on these tree pairs individually not only are the trees resolved in an incompatible manner the transitivity of the MCCs is also broken. If no reassortment has occurred between leaves `A` and `B` in `tree a` and `tree c` and no reassortment has occurred between these leaves in `tree b` and `tree c`, reassortment cannot have occurred between `A` and `B` in `tree a` and `tree b`. However, there has clearly been a reassortment event between trees `tree a` and `tree b`. 

However, fixing resolution issues as described above does not necessarily fix transitivity. `TreeKnit` uses simulated annealing and removes branches at random. This means that even if `tree c` is now resolved according to `tree a` and we infer that a reassortment event has happened between trees `tree a` and `tree b` as well as between trees `tree c` and `tree b` the MCCs that `TreeKnit` infers might be inconsistent with each other. For example look at the following MCCs:
```
{ 
    "MCC_dict" : {
        "1": { 
            "trees":["a", "b"],
            "mccs": [["A"],["B","C"]]
            },
        "2": { 
            "trees":["a", "c"],
            "mccs": [["A","B","C"]]
            },
        "3": { 
            "trees":["b", "c"],
            "mccs": [["A","B"],["C"]]
        }
    }
}
```
These MCCs show that no reassortment has occurred between leaves `B` and `C` in `tree a` and `tree b` and no reassortment has occurred between these leaves in `tree_a` and `tree_c`. Thus, for transitivity to hold reassortment cannot have occurred between `B` and `C` in `tree b` and `tree c`. But this is not the case.
Such inconsistencies make it impossible to visualize an ARG. Currently, we are not able to fully fix such incompatibilities. 

## Parallel MultiTreeKnit

For 4 or more trees the `--parallel` flag can be used to run `MultiTreeKnit` in parallel. By running trees in parallel we can improve run time from order $k^2$ to $k$, where $k$ is the number of trees. Running pair MCCs in parallel while keeping consistency conditions requires an ordering of jobs. This ordering can be easily determined from the recursive order graph. 

![plot](./Pictures/resolution_order.png) 

For example, in order to compute the MCCs of `tree 1` and `tree 4` or $T_{1,4}$ we must have calculated $T_{1,2}$ and $T_{1,3}$. In order to compute the MCCs of `tree 2` and `tree 3` we also need $T_{1,2}$ and $T_{1,3}$ but we do not need to know $T_{1,4}$. Thus, $T_{1,4}$ and $T_{2,3}$ can be calculated at the same time, $T_{2,4}$ must wait for $T_{1,4}$ and $T_{2,3}$ to finish, and $T_{3,4}$ must in turn wait for $T_{2,4}$. Thus, the arrows in the recursive order graph also determine the work flow in parallel computing. As can be seen the longest path is from $T_{1,2}$ to $T_{1,4}$ and then down to $T_{3,4}$, and is of length $2k - 2$, whereas there are a total of $k(k-1)/2$ pairs which would mean a quadratic runtime without parallelization. 

Per default the `--parallel` flag is set to false. When using `TreeKnit` on a larger computing cluster it is best practice to set the desired number of cluster nodes to 1, and cpu number on that node as high as desired. 
