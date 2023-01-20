```@meta
CurrentModule = TreeKnit
DocTestSetup  = quote
    using TreeTools
end	
```

# Functions

```@index
Pages = ["functions.md"]
```

## Function used for the CLI
```@docs
TreeKnit.treeknit
```

## Computing Maximal Compatible Clades (MCCs) for pairs
### Main functions 
```@docs
run_treeknit!
```

### For pairs of trees only
```@docs
naive_mccs
TreeKnit.runopt
```


## Resolving trees
### Using topology
```@docs
resolve!(::TreeTools.Tree, ::TreeTools.Tree)
```