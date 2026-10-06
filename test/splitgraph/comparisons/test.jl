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

@testset "Split layout" begin
	L = SplitGraph.SplitLayout(g)
	@test sort(L.bit) == 0:(n-1)
	for (a, node) in enumerate(g.internals)
		bits = [64*(w-1) + b for w in axes(L.clade, 1) for b in 0:63 if (L.clade[w, a] >> b) & 1 == 1]
		@test sort(bits) == sort(L.bit[node.conf])
		@test all(iszero, L.clade[setdiff(axes(L.clade, 1), L.wlo[a]:L.whi[a]), a])
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
	st = SplitGraph.EnergyState(SplitGraph.SplitLayout(g), ones(Bool, n), resolve)
	@test st.E == compute_energy(st.conf, g)
	nremoved = 0
	F = st.E + γ*nremoved
	for T in (1., 0.1)
		for _ in 1:300
			F, nremoved, _ = SplitGraph.mcmcstep!(st, F, nremoved, T, γ)
			@test st.E == compute_energy(st.conf, g)
			@test nremoved == count(!, st.conf)
			@test F == st.E + γ*nremoved
		end
	end
	SplitGraph.set_resolve(true)
end

@testset "Incremental energy of three trees (resolve=$resolve)" for resolve in (false, true)
	trees = [parse_newick_string(nwk) for nwk in (
		"(((A,B),(C,(D,E))),((F,G),(H,(I,J))),K,L);",
		"(((A,C),(B,(D,K))),((F,(G,L)),(H,I)),J,E);",
		"((A,(B,C,D)),(E,F,G),((H,I),(J,K,L)));",
	)]
	g3 = trees2graph(trees...)
	SplitGraph.set_resolve(resolve)
	rng = MersenneTwister(3)
	st = SplitGraph.EnergyState(SplitGraph.SplitLayout(g3), ones(Bool, length(g3.leaves)), resolve)
	for _ in 1:2000
		E, count, term = st.E, copy(st.count), copy(st.term)
		j = rand(rng, eachindex(st.conf))
		@test SplitGraph.flip!(st, j) == compute_energy(st.conf, g3)
		if rand(rng) < 0.5
			SplitGraph.undo!(st)
			@test (st.E, st.count, st.term) == (E, count, term)
			@test st.E == compute_energy(st.conf, g3)
		end
	end
	SplitGraph.set_resolve(true)
end

# Annealing at one temperature as it was written before the incremental energy: the full
# energy after every flip, and a search of the list at every step of minimal free energy.
function doMCMC_reference(g, conf, M; T, γ)
	_conf = copy(conf)
	F = compute_energy(_conf, g) + γ*count(!, _conf)
	Fmin = F
	oconf = [copy(_conf)]
	for _ in 1:M
		i = rand(1:length(_conf))
		_conf[i] = !_conf[i]
		Fnew = compute_energy(_conf, g) + γ*count(!, _conf)
		if Fnew < F || exp(-(Fnew-F)/T) > rand()
			F = Fnew
		else
			_conf[i] = !_conf[i]
		end
		if F < Fmin
			Fmin = F
			oconf[1] .= _conf
			deleteat!(oconf, 2:length(oconf))
		end
		if F == Fmin && !in(_conf, oconf)
			push!(oconf, copy(_conf))
		end
	end
	return oconf, _conf, Fmin
end

@testset "Annealing matches the full-energy reference (resolve=$resolve)" for resolve in (false, true)
	SplitGraph.set_resolve(resolve)
	rng = MersenneTwister(4)
	for (T, p) in ((1., 0.), (0.2, 0.1), (0.05, 0.3))
		conf = Vector{Bool}(rand(rng, n) .>= p)
		Random.seed!(5)
		expected = doMCMC_reference(g, conf, 400; T, γ = 2.)
		Random.seed!(5)
		@test SplitGraph.doMCMC(g, conf, 400; T, γ = 2.) == expected
	end
	SplitGraph.set_resolve(true)
end

# Likelihood of a configuration as it was written before the split layout: climbs that
# count kept leaves in each split, and branch lengths from `TreeTools.divtime`.
function conf_likelihood_reference(conf, g, seq_lengths, trees)
	L, Z = 0., 0.
	for (i, s) in enumerate(conf)
		for k1 in 1:g.K
			tn1 = trees[k1].lleaves[g.labels[i]]
			if s
				a1, ta1 = g.leaves[i].anc[k1], tn1.anc
				while SplitGraph.is_trivial_split(a1.conf, conf) && !a1.isroot
					a1, ta1 = a1.anc, ta1.anc
				end
				for k2 in (k1+1):g.K
					tn2 = trees[k2].lleaves[g.labels[i]]
					a2, ta2 = g.leaves[i].anc[k2], tn2.anc
					while SplitGraph.is_trivial_split(a2.conf, conf) && !a2.isroot
						a2, ta2 = a2.anc, ta2.anc
					end
					if SplitGraph.get_resolve() ? SplitGraph.are_equal_with_resolution(g, a1, a2, conf) :
						SplitGraph.are_equal(a1, a2, conf)
						L += SplitGraph.branch_likelihood(
							divtime(tn1, ta1), divtime(tn2, ta2), seq_lengths[k1], seq_lengths[k2]
						)
						Z += 1.
					end
				end
			else
				for k2 in 1:g.K
					tn2 = trees[k2].lleaves[g.labels[i]]
					L -= SplitGraph.branch_likelihood(
						divtime(tn1, tn1.anc), divtime(tn2, tn2.anc), seq_lengths[k1], seq_lengths[k2]
					)
					Z += 1.
				end
			end
		end
	end
	return L/max(Z, 1)
end

@testset "Likelihood matches the reference bit for bit (resolve=$resolve)" for resolve in (false, true)
	SplitGraph.set_resolve(resolve)
	rng = MersenneTwister(6)
	layout = SplitGraph.SplitLayout(g)
	for p in (0., 0.05, 0.3), _ in 1:5
		conf = Vector{Bool}(rand(rng, n) .>= p)
		expected = conf_likelihood_reference(conf, g, [1000, 1500], [t1, t2])
		@test SplitGraph.conf_likelihood(conf, g, [1000, 1500], [t1, t2]; layout) === expected
	end
	SplitGraph.set_resolve(true)
end

@testset "Distance to an ancestor equals the tree distance" begin
	for t in (t1, t2), node in values(t.lnodes), a in (node, node.anc, t.root)
		isnothing(a) && continue
		@test SplitGraph.distance_to_ancestor(node, a) === TreeTools.distance(node, a)
	end
end
