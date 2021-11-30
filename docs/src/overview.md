# Overview

## Using the CLI

`RecombTools` offers a simple CLI script: `treeknit`. 
In short, it takes two trees as input, passed as [Newick](https://en.wikipedia.org/wiki/Newick_format) files, and returns an Ancestral Reassortment Graph. 

!!! info "Compile time" 
    Julia is compiled *just in time*, meaning that functions are compiled when called for the first time inside a julia session. For this reason, some compilation will take place each time the  `treeknit` script is called, leading to a $\sim 1$ second overhead. If `Treeknit` is to be applied to many pairs of trees, it will be faster to use it from a Julia session. 

### Input

Inputs to `treeknit` are two files containing trees in Newick format. 
The only strict condition on the trees is that they share their leaf nodes. 
Example: 
```
recombtools treeknit tree1.nwk tree2.nwk
```

!!! warning "Insignificant branches"
    Tree builders sometimes introduce branches of insignificant length in order to resolve polytomies and obtain binary trees. Since `RecombTools` relies on topological differences between trees, it is important to remove these branches prior to passing the trees to `treeknit`. This can be done by using only branches with high bootstrap value (typically, $>75$), or by removing branches shorter than, *e.g.*, $(L/2)^{-1}$, where $L$ is the length of the sequences. 

!!! warning "Rooting the trees"
    The result `RecombTools` depend how the trees are rooted. It is important that the two trees are rooted in a consistent way. We recommend using the same outgroup for rooting both trees.

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

## Using from Julia
