# TreeKnit

## Installation

You can simply install *TreeKnit* from the url of the github repo using the julia package manager. 
```julia
using Pkg
Pkg.add("url="https://github.com/PierreBarrat/TreeKnit#master")
```

You should now be able to use `using TreeKnit` from inside julia. 
!!! info "TreeTools package"
    If you are going to use `TreeKnit` from inside a julia session, you will very likely need the `TreeTools` package to read Newick tree files in a format that `TreeKnit` takes as input. You can get `TreeTools` by typing Pkg.add("TreeTools") from a julia console. 

To use the CLI (Linux/Mac users), build the package by calling 
```julia
using Pkg
Pkg.build("TreeKnit")
```

This will add executable scripts to your `~/.julia/bin` folder. 
Simply add this folder to your path to call the script, *e.g.* `export PATH="$HOME/.julia/bin:$PATH"`. 
You should now be able to call, *e.g.*, `treeknit --help`

