# RecombTools

## Installation

*RecombTools* relies on the non-registered julia package *TreeTools*, that you have to install first. 
  After that, simply install *RecombTools* from the url of this github repo: 
```julia
using Pkg
Pkg.add(url="https://github.com/PierreBarrat/TreeTools#master")
Pkg.add(url="https://github.com/PierreBarrat/RecombTools#master")
```

You should now be able to use `using RecombTools`.  

To use the CLI (Linux/Mac users), build the package by calling 
```julia
Pkg.build("RecombTools")
```

This will add executable scripts to your `~/.julia/bin` folder. 
Simply add this folder to your path to call the script, *e.g.* `export PATH="$HOME/.julia/bin:$PATH"`. 
You should now be able to call, *e.g.*, `recombtools treeknit --help`