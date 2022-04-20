# MultiTreeKnit.jl
Extension of the julia package TreeKnit for inference of Ancestral Reassortment Graphs of multiple segments, the goal is to extend the package [TreeKnit](https://pierrebarrat.github.io/TreeKnit.jl), which performs inference of Ancestral Reassortment Graphs for segmented genomes (typically, human influenza) using two segment trees to a multiple segments. This is important as influenza is comprised of 8 segments and thus enabling the inference of ARGs between all segments could be useful in understanding influenza evolution processes. 

This folder just contains code for running a benchmark version of multi-tree TreeKnit. 
In this approach trees (with the same leaves) are given as input in a specific order and this is the order used to resolve polytomies and infer MCCs. The order that pairs are resolved in is shown in the picture below 

![plot](./Pictures/resolution_order.png)

Here T_{1,2} corresponds to resolving trees 1 and 2 using each other and then generating a list of MCCs using TreeKnit, these resolved trees are then used to calculate the MCCs of the next neighboring tree pairs (which are connected by arrows). This means that after calculating the MCCs for tree 1 and tree 2 and resolving their polytomies using each other, instead of using the original tree 1 the resolved tree 1 is then further used when calculating the MCCs between tree 1 and tree 3. 

The complete description of the TreeKnit algorithm is available at
> TreeKnit: Inferring Ancestral Reassortment Graphs of influenza viruses   
> Pierre Barrat-Charlaix, Timothy G. Vaughan, Richard A. Neher
> bioRxiv 2021.12.20.473456; doi: https://doi.org/10.1101/2021.12.20.473456
