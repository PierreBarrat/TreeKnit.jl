using TreeKnit

nwk1 = "../test/resolving/tree1.nwk"
nwk2 = "../test/resolving/tree2.nwk"

TreeKnit.treeknit(nwk1, nwk2)

 MCCs = [
 	["H"],
 	["A1", "A2", "B", "C", "D1", "D2", "E", "F", "G"]
 ]
