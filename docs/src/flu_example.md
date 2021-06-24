# Flu example

*RecombTools* includes the small submodule *Flu* that contains convenient functions when working with influenza trees. 
  To infer MCCs from a pair of influenza trees, we first use the `Flu.read_flu_trees` function that takes a dictionary pointing to newick files as input: 

```@example flu
using RecombTools
tree_files = Dict(
	"ha" => pathof(RecombTools) * "examples/tree_h3n2_ha.nwk",
	"na" => pathof(RecombTools) * "examples/tree_h3n2_na.nwk"
)
flutrees = Flu.read_flu_trees(tree_files)
```

`read_flu_trees` does two things: 
- read input trees and store them in a `Dict{String, Tree}`, which is the input format of `computeMCCs`
- remove "short" branches from the trees. 
  Some tree builders introduce branches of insignificant length in trees to make them binary, which conflicts with the topological approach used by `RecombTools`. 
  It is necessary to remove these branches. 
  The threshold below which a branch is removed is $\frac{1}{2L}$ where $L$ is the length of the gene sequence. 


