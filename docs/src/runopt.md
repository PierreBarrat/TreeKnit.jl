# The `runopt` function

```@setup runopt
using RecombTools
using RecombTools.SplitGraph
using TreeTools
```

- For larger trees, `opttrees` will not find all reassortments at once. After one round of `opttrees`, we're left with trees that still have incompatibilities
- Simple example with two trees. Trees not shown for space reasons, but you're encouraged to draw them. Leaves `X` and `Y` are the "obvious" reassorted ones when inspecting trees.
```@example runopt
nwk1 = "((D,(Y,((A,X),(B,C)))),E)"
nwk2 = "(Y,(E,(D,(A,(B,(C,X))))))"
t1 = node2tree(parse_newick(nwk1));
t2 = node2tree(parse_newick(nwk2));
trees = Dict(1=>t1, 2=>t2)
```