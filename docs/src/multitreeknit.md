# [MultiTreeKnit] (@id multitreeknit)

`TreeKnit` can be used to infer recombination events between multiple tree pairs. When recombination events between three or more trees should be infered `TreeKnit` uses a recursive inference strategy to return consistent results. 

It is possible to run `TreeKnit` on all tree pairs individually, however this can lead to multiple issues, for example:
- Trees could be resolved inconsistently. By default `TreeKnit` resolves trees when searching for recombination events, this is especially important for influenza where resolution is often low and not resolving trees can lead to much higher rates of rcombination being infered. However, this could also mean that when the MCCs of `tree1` and `tree2` and the MCCs of `tree1` and `tree3` are computed `tree1` might be resolved differently with `tree2` than with `tree3`. Not allowing this recombination event information to be used together.
- The final MCCs might be inconsistent with each other. Recombination is transitive. If no recombination has occured between leaves `a` and `b` in `tree1` and `tree2` and no recombination has occured between these leaves in `tree1` and `tree3`, recombination cannot have occured between `a` and `b` in `tree2` and `tree3`. However, due to the fact that `TreeKnit` uses simulated annealing and removes branches at random this might still happen (see example). Such inconsistencies make it impossible to visualize an ARG.

## Consistent Resolution
To avoid inconsistent resolution, the tree order is used to resolve tree polytomies consistently and infer MCCs. The order that pairs are resolved in is shown in the picture below 

![plot](./Pictures/resolution_order.png)

Here $T_{1,2}$ corresponds to resolving trees 1 and 2 using each other and then generating a list of MCCs using `TreeKnit`, these resolved trees are then used to calculate the MCCs of the next neighboring tree pairs (which are connected by arrows). This means that after calculating the MCCs for `tree1` and `tree2` and resolving their polytomies using each other, instead of using the original `tree1` the `resolved_tree1` is then further used when calculating the MCCs between `tree1` and `tree3`. 

Note that even if the `resolved_tree1` is further resolved with `tree3` the MCCs found between `tree1` and `tree2` will still hold for this twice resolved `tree1`.

## Consistent MCCs
To avoid inconsistent MCCs we use two approaches, one we give previously infered MCCs as a constraint to the simulated annealing. For example, assume that in `tree1` and `tree2` as well as in `tree1` and `tree3` the leaves `a` and `b` were both found to be in an MCC together and no recombination event to have occured between them. The branches between leaves `a` and `b` are marked as shared and their removal is assigned a higher cost (`oa.constraint_cost `$= 2*\gamma$) than removing a node that is not on a shared branch. In this manner we push simulated annealing to find MCCs between `tree2` and `tree3` that are consistent with the other MCCs. 

However, as simulated annealing is a stochastic process it is still possible that MCCs will be found that are inconsistent with each other. Therefore, if an ARG of all trees is desired the `--force-consist` flag can be used to make all MCCs consistent with each other. Potentially by splitting MCCs. 


