# Flu example

*TreeKnit* includes the small submodule *Flu* that contains convenient functions when working with influenza trees. 
  To infer MCCs from a pair of influenza trees, we first use the `Flu.read_flu_trees` function that takes a dictionary pointing to newick files as input: 

```@example flu
using TreeKnit
tree_files = Dict(
	"ha" => dirname(pathof(TreeKnit)) * "/../examples/tree_h3n2_ha_2013-09.nwk",
	"na" => dirname(pathof(TreeKnit)) * "/../examples/tree_h3n2_na_2013-09.nwk"
)
flutrees = Flu.read_flu_trees(tree_files)
```

`read_flu_trees` does two things: 
- read input trees and store them in a `Dict{String, Tree}`, which is the input format of `computeMCCs`
- remove "short" branches from the trees. 
  Some tree builders introduce branches of insignificant length in trees to make them binary, which conflicts with the topological approach used by `TreeKnit`. 
  It is necessary to remove these branches. 
  The threshold below which a branch is removed is $\frac{1}{2L}$ where $L$ is the length of the gene sequence. 

We can now simply infer MCCs for these trees:

```@example flu
computeMCCs(flutrees)["ha", "na"]
```

