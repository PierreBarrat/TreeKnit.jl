using Test
using TreeKnit.SRG
using TreeTools

## Simple case
nwk1 = "(A:3,(B:2,C:2)BC:1)R;"
nwk2 = "((A:3,B:3)AB:1,C:4)R;"

t1 = parse_newick_string(nwk1);
t2 = parse_newick_string(nwk2);

MCCs = [["A","B"], ["C"]]

arg, rlm, lm1, lm2 = SRG.arg_from_trees(t1, t2, MCCs)


## Case with a lot of singletons
nwk1 = "(((((A:2,X3)AX3:1,X4)AX4:2,X5)AX5:2,B)AB:0.,X1,X2):0.R;"
nwk2 = "((((A:1,X2)AX2:3,X1)AX1:3,B)AB:0.,X3,X4,X5):0.R;"

t1 = parse_newick_string(nwk1);
t2 = parse_newick_string(nwk2);

MCCs = [["A","B"],["X1"],["X2"],["X3"],["X4"],["X5"]]

arg, rlm, lm1, lm2 = SRG.arg_from_trees(t1, t2, MCCs)

## Case with two different roots
nwk1 = "(((A:1,B:1)AB:1,C:2)ABC:1,D:3)R;"
nwk2 = "((A:3,D:3)AD:1,(B:3,C:3)BC:1)R;"

t1 = parse_newick_string(nwk1);
t2 = parse_newick_string(nwk2);

MCCs = [["A","D"], ["B"], ["C"]]

arg, rlm, lm1, lm2 = SRG.arg_from_trees(t1, t2, MCCs)
