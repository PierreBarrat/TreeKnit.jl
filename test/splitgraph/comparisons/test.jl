using Random
using Test
using TreeKnit
using TreeKnit.SplitGraph
using TreeTools

t1 = read_tree("$(dirname(pathof(TreeKnit)))/..//test/NYdata/tree_ha.nwk")
t2 = read_tree("$(dirname(pathof(TreeKnit)))/..//test/NYdata/tree_na.nwk")
g = trees2graph(t1, t2)
n = length(g.leaves)

kept(a, conf) = Set(i for i in a.conf if conf[i])

@testset "Membership vectors of split nodes" begin
	for a in g.internals
		@test length(a.mask) == n
		@test findall(a.mask) == a.conf
	end
end

@testset "Split comparisons match set operations" begin
	rng = MersenneTwister(1)
	nodes = g.internals
	for p in (0., 0.2, 0.6)
		conf = rand(rng, n) .>= p
		for _ in 1:2000
			a1, a2 = rand(rng, nodes), rand(rng, nodes)
			s1, s2 = kept(a1, conf), kept(a2, conf)
			@test SplitGraph.is_contained(a1, a2, conf) == issubset(s1, s2)
			@test SplitGraph.are_equal(a1, a2, conf) == (s1 == s2)
			@test SplitGraph.are_disjoint(a1, a2, conf) == isempty(intersect(s1, s2))
		end
	end
end

@testset "Annealing step tracks energy and removed leaves (resolve=$resolve)" for resolve in (false, true)
	SplitGraph.set_resolve(resolve)
	Random.seed!(2)
	γ = 2.
	conf = ones(Bool, n)
	E = compute_energy(conf, g)
	nremoved = 0
	F = E + γ*nremoved
	for T in (1., 0.1)
		for _ in 1:300
			E, F, nremoved = SplitGraph.mcmcstep!(conf, g, E, F, nremoved, T, γ)
			@test E == compute_energy(conf, g)
			@test nremoved == count(!, conf)
			@test F == E + γ*nremoved
		end
	end
	SplitGraph.set_resolve(true)
end
