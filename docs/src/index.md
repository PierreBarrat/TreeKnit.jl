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

Should you want to simulate ARGs (*i.e.* multiple trees with reassortment events), you might find it useful to get the *ARGTools* package: 
```julia
using Pkg
Pkg.add(url="https://github.com/PierreBarrat/ARGTools#master")
```