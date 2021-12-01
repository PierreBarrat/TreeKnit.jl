@testset "Creating" begin
	sn = SplitNode()
	ln = LeafNode([], 1)
	@test true
end

t1 = read_tree("$(dirname(pathof(TreeKnit)))/..//test/splitgraph/basic/tree1.nwk")
t2 = read_tree("$(dirname(pathof(TreeKnit)))/..//test/splitgraph/basic/tree2.nwk")
g = trees2graph(t1,t2)

@testset "Basic SplitGraph construction" begin
	@testset "Leaves" for i in 1:length(leaves(t1))
		@test g.leaves[i].conf == Int[i]
		@test g.labels_to_int[g.labels[i]] == i
	end
	#
	leaf = "C"
	idx = g.labels_to_int[leaf]
	@test g.leaves[idx].anc[1].conf == sort(Int[g.labels_to_int[x] for x in ["C","D"]])
	@test g.leaves[idx].anc[1].anc.conf == sort(Int[g.labels_to_int[x] for x in ["C","D","E"]])
	@test g.leaves[idx].anc[1].anc.anc.conf == sort(Int[g.labels_to_int[x] for x in ["C","D","E","A","B"]])
	@test g.leaves[idx].anc[2].conf == sort(Int[g.labels_to_int[x] for x in ["C","B"]])
	@test g.leaves[idx].anc[2].anc.conf == sort(Int[g.labels_to_int[x] for x in ["C","B","A"]])
	@test g.leaves[idx].anc[2].anc.anc.conf == sort(Int[g.labels_to_int[x] for x in ["C","D","E","A","B"]])
end

@testset "Basic compute_energy" begin
	conf = ones(Bool, length(g.leaves))
	@test compute_energy(conf, g) == 5
	@testset for i in Iterators.filter(!=(3), 1:5)
		conf[i] = false
		@test compute_energy(conf, g) == 4
		conf[i] = true
	end
	conf = Bool[1, 1, 0, 1, 1]
	@test compute_energy(conf, g) == 0
	@testset for i in 1:5
		conf[i] = false
		@test compute_energy(conf, g) == 0
	end
end

# For exactly matching trees
t1 = node2tree(parse_newick("((A,B),(C,D))"))
t2 = node2tree(parse_newick("((A,B),(C,D))"))
g = trees2graph(t1,t2);

@testset "Compute energy for identical trees" begin
	conf = ones(Bool, length(g.leaves))
	@test compute_energy(conf, g) == 0
	@test SplitGraph.count_mismatches(t1,t2) == 0
end
