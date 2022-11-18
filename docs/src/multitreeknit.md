# [MultiTreeKnit] (@id multitreeknit)

`TreeKnit` can be used to infer recombination events between multiple tree pairs. When recombination events between three or more trees should be inferred `TreeKnit` uses a recursive inference strategy to:
- resolve each tree using information from all other trees in a consistent manner,
- identify shared regions of each tree pair (MCCs),
- find recombination events between tree pairs.
The following example illustrates the benefits of this inference strategy over running `TreeKnit` individually on all tree pairs. Note however, that MultiTreeKnit may still return MCCs that are inconsistent with each other, which prevents the construction of an ARG. Therefore, we do not reconstruct an ARG for more than two trees. 

Furthermore, MultiTreeKnit has been optimized for using multiple reassorting segments to infer parameters of pathogen evolution. We optimized for sensitivity in MCC inference and polytomy resolution, minimizing the false positive rate for shared branches and new splits, but increasing the false negative rate, i.e. typically inferring too many reassortment events. 

![plot](./Pictures/MultiTK_example.png)
## Consistent Resolution

When `TreeKnit` is run individually on all tree pairs trees could be resolved inconsistently. By default `TreeKnit` resolves trees when searching for recombination events, this is especially important for influenza where resolution is often low and not resolving trees can lead to much higher rates of recombination being inferred. 

For the example above this would mean that when the MCCs of `tree a` and `tree c` are computed no recombination events would be found as `tree c` would be resolved according to `tree a` (i.e. a branch would be introduced above leaves `A` and `B`). However, when the MCCs of `tree b` and `tree c` are computed `tree c` would be resolved according to `tree b` and a branch would be introduced above leaves `B` and `C`. These two different resolutions of `tree c` are incompatible with each other and do not allow us to use this recombination event information together.

To avoid inconsistent resolution, the tree order is used to resolve tree polytomies consistently and infer MCCs. The order that pairs are resolved in is shown in the picture below 

![plot](./Pictures/resolution_order.png)

Here $T_{1,2}$ corresponds to resolving `tree 1` and `tree 2` using each other and then generating a list of MCCs using standard `TreeKnit`. These resolved trees are then used to calculate the MCCs of the next neighboring tree pairs (which are connected by arrows). This means that after calculating the MCCs for `tree 1` and `tree 2` and resolving their polytomies using each other, instead of using the original `tree 1` the `resolved tree 1` is then further used when calculating the MCCs between `tree 1` and `tree 3`. Leading to a consistent tree resolution.

Furthermore, we use strict resolve when resolving trees with MultiTreeKnit. Strict resolve will only resolve polytomies if the location of each branch in the new split can be fully determined by the other tree. Using strict resolve prevents the introduction of incorrect splits into the trees, this is especially important when resolved trees are used downstream for inference as these splits could prevent the simulated annealing from converging. 

At the end of the sequential inference on all tree pairs, each tree will be resolved as much as possible using each other tree. However, the output MCCs might not be consistent with each other and might not fulfill the necessary transitivity requirements to create an ARG. To reduce the number of such potential incompatibilities we run a final round where we re-infer all tree pair MCCs without resolving trees this will prevent topological incompatibilities and ensure that all inferred MCCs are compatible with the output resolved trees, however it may not fix incompatibilities that are due to transitivity, therefore we do not produce an ARG for more than two trees. 

## Consistent MCCs

When `TreeKnit` is run individually on all tree pairs the final MCCs might be inconsistent with each other. 
As can be seen in the example when `TreeKnit` is run on these tree pairs individually not only are the trees resolved in an incompatible manner the transitivity of the MCCs is also broken. If no recombination has occurred between leaves `A` and `B` in `tree a` and `tree c` and no recombination has occurred between these leaves in `tree b` and `tree c`, recombination cannot have occurred between `A` and `B` in `tree a` and `tree b`. However, there has clearly been a recombination event between trees `tree a` and `tree b`. 

However, fixing resolution issues as described above does not necessarily fix transitivity. `TreeKnit` uses simulated annealing and removes branches at random. This means that even if `tree c` is now resolved according to `tree a` and we infer that a recombination event has happened between trees `tree a` and `tree b` as well as between trees `tree c` and `tree b` the MCCs that `TreeKnit` infers might be inconsistent with each other. For example look at the following MCCs:
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
These MCCs show that no recombination has occurred between leaves `B` and `C` in `tree a` and `tree b` and no recombination has occurred between these leaves in `tree_a` and `tree_c`. Thus, for transitivity to hold recombination cannot have occurred between `B` and `C` in `tree b` and `tree c`. But this is not the case.
Such inconsistencies make it impossible to visualize an ARG.

To decrease the number of inconsistencies in `MultiTreeKnit` we give previously inferred MCCs as a constraint to the simulated annealing. For example, assume that in `tree a` and `tree b` as well as in `tree a` and `tree c` the leaves `B` and `C` were both found to be in an MCC together and no recombination event to have occurred between them. The branches between leaves `B` and `C` are marked as shared in trees `tree b` and `tree c` and their removal is assigned a higher cost (`oa.constraint_cost `$= 2*\gamma$) than removing a node that is not on a shared branch. In this manner we push the simulated annealing to find MCCs between `tree b` and `tree c` that are consistent with the other MCCs, however we still enable it to split such branches if it leads to a more optimal energy cost if, for example, the previously inferred MCCs were inaccurate. By default we run two rounds of `MultiTreeKnit`. In the first round the initial tree pairs have no constraints but in the second round they also have constraints from other tree pairs, enabling optimal information transfer. In the final round we additionally do not resolve polytomies any further, as we assume by now each tree has been resolved as much as possible using every other tree. This prevents topological incompatibilities and ensures that all inferred MCCs are compatible with the output resolved trees.

However, even using these approaches MCCs may still be inconsistent as simulated annealing is a stochastic process. Inconsistent MCCs cannot be viewed together in a ARG. Currently, we are not able to fully fix such incompatibilities. 

## Parallel MultiTreeKnit

For 4 or more trees the `--parallel` flag can be used to run `MultiTreeKnit` in parallel. By running trees in parallel we can improve run time from order $k^2$ to $k$, where $k$ is the number of trees. Running pair MCCs in parallel while keeping consistency conditions requires an ordering of jobs. This ordering can be easily determined from the recursive order graph. 

![plot](./Pictures/resolution_order.png) 

For example, in order to compute the MCCs of `tree 1` and `tree 4` or $T_{1,4}$ we must have calculated $T_{1,2}$ and $T_{1,3}$. In order to compute the MCCs of `tree 2` and `tree 3` we also need $T_{1,2}$ and $T_{1,3}$ but we do not need to know $T_{1,4}$. Thus, $T_{1,4}$ and $T_{2,3}$ can be calculated at the same time, $T_{2,4}$ must wait for $T_{1,4}$ and $T_{2,3}$ to finish, and $T_{3,4}$ must in turn wait for $T_{2,4}$. Thus, the arrows in the recursive order graph also determine the work flow in parallel computing. As can be seen the longest path is from $T_{1,2}$ to $T_{1,4}$ and then down to $T_{3,4}$, and is of length $2k - 2$, whereas there are a total of $k(k-1)/2$ pairs which would mean a quadratic runtime without parallelization. 

Per default the `--parallel` flag is set to false. When using `TreeKnit` on a larger computing cluster it is best practice to set the desired number of cluster nodes to 1, and cpu number on that node as high as desired. 
