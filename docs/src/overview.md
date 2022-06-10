# Overview

## Using the CLI

`TreeKnit` offers a simple CLI script: `treeknit`. 
In short, it takes two trees as input, passed as [Newick](https://en.wikipedia.org/wiki/Newick_format) files, and returns an Ancestral Reassortment Graph. 

!!! info "Compile time" 
    Julia is compiled *just in time*, meaning that functions are compiled when called for the first time inside a julia session. For this reason, some compilation will take place each time the  `treeknit` script is called, leading to an overhead of a few seconds. If `Treeknit` is to be applied to many pairs of trees, it will be faster to use it from a Julia session. 

### Input

Inputs to `treeknit` are two files containing trees in Newick format. 
The only strict condition on the trees is that they share their leaf nodes. 
Example: 
```
treeknit tree1.nwk tree2.nwk
```

There are other important requirements for input trees, highlighted below. 
They are not strict conditions in the sense that the algorithm will not fail if they are not met. 
However, not meeting them might result in irrelevant or meaningless output. 

!!! warning "Insignificant branches"
    Tree builders sometimes introduce branches of insignificant length in order to resolve polytomies and obtain binary trees. Since `TreeKnit` relies on topological differences between trees, and all internal nodes of input trees are interpreted as hard topological constraints. It is therefore important to remove low support branches/internal nodes prior to passing the trees to `treeknit`. This is most easily done by removing every branch that is not supported by at least one mutation, *e.g.* branches shorter than $(L/2)^{-1}$, where $L$ is the length of the sequences. Additionally, internal nodes with a low bootstrap value (typically $<75$) can be remove for additional robustness. 

!!! warning "Rooting the trees"
    The result `TreeKnit` depend how the trees are rooted. It is important that the two trees are rooted in a consistent way. We recommend using the same outgroup for rooting both trees.

### Output

Output of the inference is written a directory called `treeknit_results`. This can be changed using the `--outdir` option.   
The directory will contain:   
- the ARG, written as an extended [Newick string](https://doi.org/10.1186/1471-2105-9-532).   
- the MCCs, *i.e.* shared regions of the trees, indicated by the leaves they contain.  
- resolved trees, where polytomies have been reduced as much as possible using the knowledge of the MCCs.   
- a table with the correspondence between internal nodes of the ARG and the trees. Note that this refers to node labels of the resolved trees, which may not be the same as the ones given as input.   

### Options

The main options that you can play with are:  
- the parsimony parameter $\gamma$, `--gamma` or `-g`.   
- naive inference `--naive`. Using this flag is equivalent to setting $\gamma \rightarrow \infty$.  
- Length of sequences used to infer trees: `--seq-lengths`. These are used for likelihood test to break degeneracy between topologically equivalent MCCs.  

More details in the [options section](@ref options).

## Using from a Julia session

If `TreeKnit` has to be used on several pairs of trees and speed is important, then you should call it from a julia session directly. 
Let's see how one does this using the example directory, which contains two Newick files `tree_h3n2_ha.nwk` and `tree_h3n2_na.nwk`. 
First, read the trees: 
```@example usage_from_julia
using TreeTools
using TreeKnit
t_ha = read_tree(dirname(pathof(TreeKnit)) * "/../examples/tree_h3n2_ha.nwk")
t_na = read_tree(dirname(pathof(TreeKnit)) * "/../examples/tree_h3n2_na.nwk")
```

We now proceed in three steps: 
1. Compute MCCs for these two trees. See the [options](@ref options) or [MCCs](@ref MCCs) for more details.
2. Resolve the trees using these MCCs. 
3. Compute the ARG from the resolved trees and the MCCs. 

```@repl usage_from_julia
MCCs = computeMCCs(t_ha, t_na) # compute MCCs
rS = resolve!(t_ha, t_na, MCCs); # resolve. Output `rS` contains the resolved splits
arg, rlm, lm1, lm2 = SRG.arg_from_trees(t_ha, t_na, MCCs); # compute the ARG and mappings from tree to ARG internal nodes. 
```

To write these results to files, one could do 
```
write("arg.nwk", arg) # write the ARG
write_mccs("mcc.dat", MCCs) # write the MCCs
write_newick("tree_ha_resolved.nwk", t_ha) # write resolved HA tree
```

Note that the `treeknit` function in `src/cli.jl` follows exactly these steps. Have a look at it for a more detailed example of how to use this package from inside julia. 




