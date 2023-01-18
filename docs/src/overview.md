# Overview

## Using the CLI

`TreeKnit` offers a simple CLI script: `treeknit`. 
It takes two or more trees as input, passed as [Newick](https://en.wikipedia.org/wiki/Newick_format) files, and infers reassortment events between all *pairs* of trees. 
It does this by finding shared regions of two trees, so called [maximally compatible clades (MCCs)](@ref MCCs), within which we assume no reassortment has occurred. 
Each MCC corresponds to one reassortment event between the two trees. 

When exactly two trees are given as input, TreeKnit additionally returns an Ancestral Reassortment Graph (ARG) that combines the two trees. 
This is currently not possible for multiple trees, more information on how `TreeKnit` runs on more than two trees can be found in the [MultiTreeKnit section](@ref multitreeknit)).

!!! info "Compile time" 
    Julia is compiled *just in time*, meaning that functions are compiled when called for the first time inside a julia session. For this reason, some compilation will take place each time the  `treeknit` script is called, leading to an overhead of a few seconds. If `Treeknit` is to be applied to many pairs of trees, it will be faster to use it from a Julia session. 

### Input

Inputs to `treeknit` are files containing trees in Newick format. 
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
    `TreeKnit` depends on how the trees are rooted. It is important that the two trees are rooted in a consistent way. We recommend using the same outgroup for rooting both trees.

### Output

Output of the inference is written to the directory `treeknit_results`. This can be changed using the `--outdir` option.   
The directory will contain:    
- the MCCs, *i.e.* shared regions for paris of trees. The MCCs of all tree pairs are written in JSON format. 
  The example below shows possible MCCs for three trees of 10 leaves named `a`, `b` and `c`. 
```
{ 
    "MCC_dict" : {
        "1": { 
            "trees":["a", "b"],
            "mccs": [["5_0", "6_0"],["10_0", "1_0", "2_0", "3_0", "4_0", "7_0", "8_0", "9_0"]]
            },
        "2": { 
            "trees":["a", "c"],
            "mccs": [["1_0", "8_0", "9_0"],["10_0", "2_0", "3_0", "4_0", "5_0", "6_0", "7_0"]]
            },
        "3": { 
            "trees":["b", "c"],
            "mccs": [["5_0", "6_0"],["1_0", "8_0", "9_0"],["10_0", "2_0", "3_0", "4_0", "7_0"]]
        }
    }
}
```
- resolved trees: polytomies in each tree can be resolved using information coming from the other trees and from the MCCs (see the [resolving section](@ref resolving))
  They are written to files with the label `_resolved` added behind the original tree file name, e.g. `tree_ha_resolved.nwk`. 
- a file with all parameters used by TreeKnit during inference titled `parameters.json`.

In addition, when TreeKnit is called on exactly two trees the directory will contain a separate ARG folder with:
- the ARG, written as an extended [Newick string](https://doi.org/10.1186/1471-2105-9-532).    
- liberally resolved trees needed for the construction of an ARG, see [resolving section](@ref resolve_strict_vs_liberal). The output will be written to a file with the label `_liberal_resolved` added behind the tree label, e.g. `tree_ha_liberal_resolved.nwk`. 
- a table with the correspondence between internal nodes of the ARG and the trees. Note that this refers to node labels of the liberally resolved trees, which may not be the same as the ones given as input.   

### Options

Options that you can play with are:  
- `--gamma` or `-g`: the parsimony parameter $\gamma$ (see [here](@ref gamma)).
- `--naive`: naive inference. Using this flag is equivalent to setting $\gamma \rightarrow \infty$.  
- `--seq-lengths`: length of sequences used to infer trees. These are used for likelihood test to break degeneracy between topologically equivalent MCCs.  
- `--no-likelihood`: do not use branch length to sort the likelihood of different MCCs.
- `--no-resolve`: do not attempt to resolve trees before inferring MCCs.
- `--parallel`: run sequential MultiTreeKnit with parallelization (only used for 4 or more trees).
- `--auspice-view`: will create files that can be used to view a tanglegram of the two trees with colored maximally compatible clades in [auspice](https://docs.nextstrain.org/projects/auspice/en/stable/advanced-functionality/second-trees.html). For more information see [Visualization of MCCs in a tanglegram](@ref view_auspice).


Additional options can be used for "advanced users":
- `rounds`: number of rounds of TreeKnit. Trees are resolved at each TreeKnit run, meaning subsequent rounds can give different results. This is especially useful for more than two trees, see [MultiTreeKnit](@ref multitreeknit).
- `--no-pre-resolve`: do not compatibly resolve all trees with each other before inferring MCCs (default is to pre-resolve).
- `--resolve-all-rounds`: resolve trees before inferring pairwise MCCs in all rounds Overrides `--no-resolve`. Default for 2 trees. For more than 2 trees the default is to not resolve in the final round.
- `--liberal-resolve`: use liberal resolving method, see [here](@ref resolve_s trict_vs_liberal). In short, with this option, TreeKnit will try to resolve trees as much as possible even in ambiguous situations. It results in more resolved trees but also more wrong splits. 

More details on the options can be found throughout the documentation. 
In particular, the parsimony parameter is explained [here](@ref gamma) and resolution of trees [here](@ref resolving).

## Using from a Julia session

!!! warning "For two trees only"
    For more than two trees, the process is different. A plan for the future is to "uniformize" the way TreeKnit is called for 2 or more trees. So what's below may be deprecated then. 

If `TreeKnit` has to be used on several datasets and speed is important, then you should call it from a julia session directly. 
Let's see how one does this for a simple two tree example using the `example` directory, which contains two Newick files `tree_h3n2_ha.nwk` and `tree_h3n2_na.nwk`. 
First, read the trees: 
```@example usage_from_julia
using TreeTools
using TreeKnit
t_ha = read_tree(dirname(pathof(TreeKnit)) * "/../examples/tree_h3n2_ha.nwk", label="ha")
t_na = read_tree(dirname(pathof(TreeKnit)) * "/../examples/tree_h3n2_na.nwk", label="na")
```
!!! info "Tree Labels" 
    When computing the MCCs for tree pairs `TreeKnit` uses a tree's `label` as a unique identifier. These labels are also written to the output `.json` file. When `TreeKnit`is used via the command line the filename and/or folder name is used as a label for the input trees if these are unique, otherwise an error is thrown. When using `TreeKnit` from a julia session, trees will be assigned a random unique string as an identifying label. One can change this by calling `label!(tree, some_label)`, see the [TreeTools documentation](https://pierrebarrat.github.io/TreeTools.jl/dev/) 

We now proceed in the following steps:
1. Compute the MCCs of the two trees. By default, if `TreeKnit` receives two trees as input, these will be resolved with each other prior to inference unless `--no-resolve` is activated. See the [options](@ref options) or [MCCs](@ref MCCs) for more details, additionally see [MultiTreeKnit](@ref multitreeknit) to learn how this is handled when more than two trees are given as input.
2. Resolve the trees using these MCCs. 
3. Optionally the first input tree can be ladderized, and the polytomies of two trees can be sorted w.r.t to the MCCs, allowing for clearer visualization (e.g. if the tree pair should be later visualized in a dendrogram).
4. Compute the ARG from the resolved trees and the MCCs. 

```@repl usage_from_julia
MCCs = computeMCCs(t_ha, t_na) # compute MCCs
rS = resolve!(t_ha, t_na, MCCs); # resolve. Output `rS` contains the resolved splits
TreeTools.ladderize!(t_ha) #ladderize the first tree
TreeKnit.sort_polytomies!(tree1, tree2, MCC) #sort the polytomies 
arg, rlm, lm1, lm2 = SRG.arg_from_trees(t_ha, t_na, MCCs); # compute the ARG and mappings from tree to ARG internal nodes. 
```

To write these results to files, one could do 
```
write("arg.nwk", arg) # write the ARG
write_mccs("MCCs.json", MCCs) # write the MCCs
write_newick("tree_ha_resolved.nwk", t_ha) # write resolved HA tree
```

Note that the `treeknit` function in `src/cli.jl` follows exactly these steps for two trees, for multiple trees see the [MultiTreeKnit section](@ref multitreeknit). Have a look at it for a more detailed example of how to use this package from inside julia. 
