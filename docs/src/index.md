# TreeKnit

## Installation

*TreeKnit* relies on the non-registered julia package *TreeTools*, that you have to install first. 
  After that, simply install *TreeKnit* from the url of this github repo: 
```julia
using Pkg
Pkg.add(url="https://github.com/PierreBarrat/TreeTools#master")
Pkg.add(url="https://github.com/PierreBarrat/TreeKnit#master")
```

You should now be able to use `using TreeKnit`.  

To use the CLI (Linux/Mac users), build the package by calling 
```julia
Pkg.build("TreeKnit")
```

This will add executable scripts to your `~/.julia/bin` folder. 
Simply add this folder to your path to call the script, *e.g.* `export PATH="$HOME/.julia/bin:$PATH"`. 
You should now be able to call, *e.g.*, `recombtools treeknit --help`