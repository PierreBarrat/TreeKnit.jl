# [Options](@id options)

`run_treeknit!` uses a topogy based heuristic optimization to find maximally compatible clades. 
  Options to this heuristic are provided through the `OptArgs` object, that is optionally passed as a second argument.
  Essential options are detailed here. 	

Note that when the `OptArgs()` object is initialized it's default parameters are the same as those for running `TreeKnit` for 2 trees. When instead the constructor `OptArgs(K::Int)` is called, the default parameters for running `TreeKnit` with `K` trees are returned. 

The default parameters for `K=2` trees are the parameters for the `:better_MCCs` method, whereas the default parameters for `K>2` trees are the parameters for the `:better_trees` method. See [TreeKnit Methods](@ref TK_method_options).

## [Parsimony parameter](@id gamma)
The heuristic method used by *TreeKnit* tries to prune consistent clades from a pair of trees in order to increase a compatibility score between other clades. 
  Pruning a clade is interpreted as fixing a reassortment right above it, while increasing the compatibility between remaining clades removes reassortments. 
  A purely parsimonious heuristic should thus give the same weight to fixing a reassortment through pruning a clade and fixing one incompatibility in the trees. 

Here, we assign the score $\gamma$ to each pruned clade, and count as $1$ each incompatibility fixed for the remaining clades. 
  For $\gamma=1$, we obtain a parsimonious method that attempts to minimize the overall number of reassortments. 
  For higher values, pruning a clade must "fix" at least $\gamma$ incompatibilities to be considered a good move, making the obtained MCCs less parsimonious. 
  For infinite $\gamma$, pruning clades is impossible, and we fall back on the [naive estimation of MCCs](@ref naive_mccs). 
  
The example below illustrates the difference between different $\gamma$ values: 
```@example gamma1
using TreeKnit # hide
t1 = parse_newick_string("((((A,B),C),D),E);")
t2 = parse_newick_string("((((D,B),E),A),C);") # Same topology, but shuffled leave
nothing # hide
```
Here, pruning the two leaves `(A,C)` or `(D,E)` results in compatible trees (resp. `((B,D),E)` and `((A,B),C)`). 
  These moves each have a cost 2$\gamma$ (removing 2 clades), but bring us from trees with 5 incompatibilities to 0. 
  They will only be accepted if $\gamma \leq 2.5$. 

```@repl gamma1
run_treeknit!(t1, t2, OptArgs(γ=2))
run_treeknit!(t1, t2, OptArgs(γ=3))
```

## Resolving trees with polytomies

See [Resolving](@ref resolving) for more information about tree resolution. Note we have multiple options for resolving trees:
- `pre_resolve:` if this option is passed to `run_treeknit!` (through `OptArgs`), `resolve!` will be called on all trees prior to MCC inference. 
- `resolve:` under this option `resolve!` is called on each **pair** of trees before each iteration of the MCC inference procedure. For example if treeknit is run on trees `t1`, `t2` and `t3` and the `resolve` is set to `true`, `resolve!` will be called on each combination of tree pairs (see [MultiTreeKnit](@ref multitreeknit) for more details).  Furthermore, this option allows for tree resolution during MCC inference on tree pairs and will resolve trees after each pair-wise iteration of `TreeKnit` using the inferred MCCs.

!!! info "pre-resolve and resolve" for tree pairs
    Note that`resolve=true` will also resolve trees with each other prior to inference, and thus `pre_resolve` is not needed for 2 trees if `resolve=true`, but we distinguish between the two options for consistency with $>2$ trees).

When resolving trees using MCCs we distinguish between two options (see [strict vs liberal resolution](@ref resolve_strict_vs_liberal) for a detailed description of the two methods):
- `strict=true:` only add unambiguous splits from one tree into another.
- `strict=false` or `liberal:` resolve will resolve trees as much as possible, fully resolving shared regions of the two trees, but arbitrarily choosing the location of ambiguous splits, leading to potentially wrong splits in the output trees. 

## [Degeneracy: sorting with likelihood](@id likelihood)
When several MCC decompositions are possible, degeneracy is removed by using the `likelihood_sort` option (activated by default). 
In the example below, there are three equivalent decompositions if only topology is considered: 
```@example degeneracy
using TreeKnit # hide
t1 = parse_newick_string("((A:2,B:2):2,C:4);")
t2 = parse_newick_string("(A:2,(B:1,C:1):1);")
oa = OptArgs(likelihood_sort = false)
unique([run_treeknit!(t1, t2, oa) for rep in 1:50]) # Repeating computation many times 
```

When taking branch lengths into account, this degeneracy vanishes: 
```@example degeneracy
oa = OptArgs(likelihood_sort = true)
unique([run_treeknit!(t1, t2, oa) for rep in 1:50])
```

## Verbosity

The TreeKnit cli has four level of verbosity that can be set with the `--verbosity-level <value>` option: 
- `-1` means no visible output
- `0` means relatively little output: input file names, number of MCCs found, etc ... Less than 10 lines. (Default)
- `1` gives a bit more details about what's going on: the MCC inference process for each pair of trees is detailed
- `2` gives a lot of information and is only useful for debugging and when using small trees 

The `--verbose` or `-v` flag set the verbosity to `1` instead of `0`. 

When using from a julia session, it can be useful to also set the verbosity. 
Let's take the simple case where one wants to use `run_treeknit!(tree1, tree2, OptArgs())` with some pre-loaded trees and default arguments. 
If run like this, no info will be shown. 
To trigger verbosity, create a custom logger in the following way: 

```@example verbosity
using TreeKnit, TreeTools
t1 = parse_newick_string("((A:2,B:2):2,C:4);")
t2 = parse_newick_string("(A:2,(B:1,C:1):1);")
using Logging
verbosity = 1
logger = ConsoleLogger(LogLevel(-verbosity)) # corresponds to `--verbosity-level 1` from the CLI
output = with_logger(logger) do 
  run_treeknit!(t1, t2, OptArgs())
end
```

The above works because during inference, TreeKnit emits log messages in the form `@logmsg LogLevel(0/-1/-2) msg`. 
The default Julia logger will only show messages with a positive level, and this can be set as described above. 
See the [Logging](https://docs.julialang.org/en/v1/stdlib/Logging) and the [LoggingExtras](https://github.com/JuliaLogging/LoggingExtras.jl) package for more information. 
