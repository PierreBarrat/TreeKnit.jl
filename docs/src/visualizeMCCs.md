## [Visualization of MCCs] (@id visualization)

# ARGPlots
Functionality to plot ancestral recombination graphs using MCCs and resolved trees generated in TreeKnit.

Assume the phylogenetic trees of multiple segmented genomes have been infered (these need the same terminal nodes) and all tree have been resolved in a compatible manner. A main tree can be chosen to view the ancestral history and recombination events of all other segments using the output generated by [TreeKnit](https://pierrebarrat.github.io/TreeKnit.jl). Given a list of resolved trees and their MCCs with the main tree (see output of `TreeKnit`) this function will plot the recombination sites as seen from the first tree. For example given two trees: 

![plot](./Pictures/tree1.png) 
![plot](./Pictures/tree2.png) 

The output ARG, as seen from tree1 is

![plot](./Pictures/arg.png)

The dashed line between recombination sites can be removed by setting the `draw_connections` argument to `False`. 