nwk1 = "((A,B),(D,(E,(F,C))))"
nwk2 = "((A,(B,C)),(D,E,F))"

t1 = node2tree(parse_newick(nwk1))
t2 = node2tree(parse_newick(nwk2))

g = SplitGraph.trees2graph((t1,t2))

@testset "Test are_equal_with_resolution" begin
	conf1 = ones(Bool, length(g.leaves))
	conf2 = copy(conf1)
	conf2[g.labels_to_int["C"]] = false
	#
	SplitGraph.set_resolve(false)
	@test compute_energy(conf1, g) == 6
	@test compute_energy(conf2, g) == 2
	#
	SplitGraph.set_resolve(true)
	@test compute_energy(conf1, g) == 6
	@test compute_energy(conf2, g) == 0
end

