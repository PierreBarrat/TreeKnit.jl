export compute_energy

let n::Int64=0
	global increment_n() = (n+=1)
	global reset_n() = (n=0)
	global get_n() = n
end
let resolve::Bool = true
	global set_resolve(r) = (resolve = r)
	global get_resolve() = resolve
end

"""
"""
function compute_energy(conf::Array{Bool,1}, g::Graph)
	length(conf) != length(g.leaves) && error("`conf` and `g` do not have the same length")
	E = 0
	if length(g.leaves) + length(g.internals) == 1
		return E
	end

	@inbounds for (i,s) in enumerate(conf)
		if s
			for k1 in 1:g.K
				# Ancestors in tree k1
				a1 = g.leaves[i].anc[k1]
				# If ancestor is identical to leaf for given configuration, go up
				# ie go up to the first non trivial split
				while is_trivial_split(a1.conf, conf) && !a1.isroot
					a1 = a1.anc::SplitNode
				end
				for k2 in (k1+1):g.K
					# Same for
					a2 = g.leaves[i].anc[k2]
					while is_trivial_split(a2.conf, conf) && !a2.isroot
						a2 = a2.anc
					end
					if get_resolve() && !are_equal_with_resolution(g, a1, a2, conf)
						E += 1
					elseif !get_resolve() && !are_equal(a1.conf, a2.conf, conf)
						E += 1
					end
				end
			end
		end
	end
	return E
end

"""
	 are_equal(nsplit1, nsplit2, conf)

Are internal configurations `nsplit1` and `nsplit2` equal for leaf configuration `conf`?
"""
function are_equal(nsplit1, nsplit2, conf)
	_in = max(length(nsplit1), length(nsplit2)) > 25 ? insorted : in

	@inbounds for n1 in nsplit1
		if conf[n1] && !_in(n1, nsplit2)
			return false
		end
	end

	@inbounds for n2 in nsplit2
		if conf[n2] && !_in(n2, nsplit1)
			return false
		end
	end

	return true
end



"""
	are_equal_with_resolution(g, a1, a2, conf)

From a given leaf, we reached nodes `a1` and `a2` for the two trees, and call `asplit1` and
`asplit2` the splits of these nodes.  Three cases:
- `asplit1 == asplit2`
- `asplit1 ⊂ asplit2` AND we can resolve `asplit2` to obtain `asplit1`
- `asplit2 ⊂ asplit1` AND we can resolve `asplit1` to obtain `asplit2`

How do we proceed in the second case?
To resolve, `asplit1` needs to be compatible with all splits below `asplit2`. We know it's
included in `asplit2`. We need that for every split `S` below `asplit2`, `S ⊂ asplit1` or
`S ∩ asplit1` is empty.
Note that it is *impossible* to have `asplit1 ⫅ S`, since in these case
we would have reached `S` instead of `asplit2` when going up from the leaf.

### Example

`asplit1 = ABC`, `asplit2 = ABCXYZ`.

1. If the children of `asplit2` are *e.g.* `AB`, `C`, `X`, `YZ`, then we can clearly resolve.
2. If ... are `ABX`, `C`, `YZ`, then we can't resolve because of the first split.
"""
function are_equal_with_resolution(g, a1::SplitNode, a2::SplitNode, conf)
	if are_equal(a1.conf, a2.conf, conf)
		return true
	elseif is_contained(a1.conf, a2.conf, conf) && _are_equal_with_resolution(g, a1.conf, a2, conf)
		return true
	elseif is_contained(a2.conf, a1.conf, conf) && _are_equal_with_resolution(g, a2.conf, a1, conf)
		return true
	end
	return false
end

"""
	are_equal_with_resolution(g, asplit1, asplit2, k2, conf)

For every leaf of `a2`, check that the corresponding split `S` is either disjoint from
`asplit1` or included in `asplit1`.
"""
function _are_equal_with_resolution(g::SplitGraph.Graph, asplit1, a2::SplitNode, conf)
	for c in a2.child
		S = c.conf
		if typeof(c) == SplitNode && !is_contained(S, asplit1, conf) && !are_disjoint(S, asplit1, conf)
			return false
		end
	end
	return true
end

function is_trivial_split(nodeconf, conf)
	n = 0
	@inbounds for i in nodeconf
		if conf[i]
			n += 1
		end
	end
	return n < 2
end

#=
NOTE

is_contained and are_disjoint could be made much faster with dictionaries for large arrays.
=#

"""
	is_contained(nsplit1, nsplit2, conf)

Is `nsplit1` in `nsplit2` for leaf configuration `conf`?
"""
function is_contained(nsplit1, nsplit2, conf)
	_in = max(length(nsplit2)) > 25 ? insorted : in
	@inbounds for i in nsplit1
		if conf[i] && !_in(i, nsplit2)
			return false
		end
	end
	return true
end

"""
	are_disjoint(S1, S2, conf)

Are `S1` and `S2` disjoint for configuration `conf`?
"""
function are_disjoint(S1, S2, conf)
	Ssmall, Sbig = if length(S1) > length(S2)
		S2, S1
	else
		S1, S2
	end

	_in = length(Sbig) > 25 ? insorted : in
	@inbounds for i in Ssmall
		if conf[i] && _in(i, Sbig)
			return false
		end
	end
	return true
end


"""
"""
function compute_F(conf::Array{Bool,1}, g::Graph, γ::Real; mask=Set(), constraint_cost=γ)
	E = compute_energy(conf, g)
	c = sum([Int(i ∈ mask && !conf[i]) for i in 1:length(conf)])*constraint_cost
	return E + γ*(length(conf) - sum(conf)) + c
end

"""
"""
function doMCMC(
	g::Graph, conf::Array{Bool,1}, M::Int64;
	T=1, γ=1, mask=Set(), constraint_cost=γ
)
	_conf = copy(conf)
	E = compute_energy(_conf, g)
	c = sum([Int(i ∈ mask && !conf[i]) for i in 1:length(conf)])*constraint_cost
	F = E + γ*(length(conf) - sum(conf)) + c
	Fmin = F

	oconf = [copy(_conf)]
	for m in 1:M
		E, F = mcmcstep!(_conf, g, F, T, γ; mask, constraint_cost)
		# If new minimum is found
		if F < Fmin
			Fmin = F
			oconf[1] .= _conf
			deleteat!(oconf, 2:length(oconf))
		end
		# If equal minimum is found
		if F == Fmin && !in(_conf, oconf) #mapreduce(x->x!=_conf, *, oconf)
			push!(oconf, copy(_conf))
		end
	end
	return oconf, _conf, Fmin
end

"""
	mcmcstep!(conf, g, F, T, γ; mask=Set(), constraint_cost=γ)

Perform an mcmc step by removing a node from `conf` at random.
If `mask` is specified masked nodes will be removed at cost `constraint_cost`,
per default this is set to the cost of a reassortment event.
"""
function mcmcstep!(conf, g, F, T, γ; mask=Set(), constraint_cost=γ)
	s = Set(1:length(conf))
	i = rand(s)
	conf[i] = !conf[i]
	Enew = compute_energy(conf, g)
	c = sum([Int(i ∈ mask && !conf[i]) for i in 1:length(conf)])*constraint_cost
	Fnew = Enew + γ*(length(conf) - sum(conf)) + c
	if Fnew < F || exp(-(Fnew-F)/T) > rand()
		return Enew, Fnew
	else
		conf[i] = !conf[i]
		return (round(Int64, F-γ*(length(conf) - sum(conf))-c), F)
	end
end

"""
	sa_opt(g::Graph; Trange=1.:-0.01:0.1, γ=1.05, M=1000, rep=1, resolve=true)

Call `_sa_opt` repeatedly to find a set of optimal confs.
"""
function sa_opt(
	g::Graph;
	Trange=OptArgs().Trange, γ=2., M=10, rep=1, resolve=true, mask=Set(), constraint_cost=γ
)
	set_resolve(resolve)
	#
	oconf = Any[]
	F = Int64[]
	Fmin = Inf
	nfound = 0
	for r in 1:rep
		oconf_, F_ = _sa_opt(g, γ, Trange, M; mask, constraint_cost)
		Fm = minimum(F_)
		if Fm == Fmin
			append!(oconf, oconf_)
			append!(F, F_)
			nfound += 1
		elseif Fm < Fmin
			Fmin = Fm
			oconf = oconf_
			F = F_
			nfound = 1
		end
	end
	return unique(oconf), F, nfound
end

function _sa_opt(g::Graph, γ, Trange, M; mask=Set(), constraint_cost= γ)
	reset_chance = 0.
	conf = ones(Bool, length(g.leaves))
	oconf = [copy(conf)]
	F = Float64[Inf]
	Fmin = F[1]
	p = Progress(
		length(Trange);
		dt=1.,
		desc="Simulated annealing: ",
		enabled=v(),
		barglyphs=BarGlyphs("[=> ]"),
		showspeed=true,
	)
	for T in Trange
		if rand() < reset_chance
			tmp_oconf, conf, fmin = doMCMC(g, oconf[rand(1:length(oconf))]; M, T, γ, mask, constraint_cost)
		else
			tmp_oconf, conf, fmin = doMCMC(g, conf, M; T,γ, mask, constraint_cost)
		end
		append!(F,fmin)
		# If a better conf is found than all configurations in oconf
		# (which is the min of `f` from doMCMC), completely replace oconf
		if fmin < Fmin
			oconf = tmp_oconf
			Fmin = fmin
		# If equally good confs have been found
		elseif fmin == Fmin
			append!(oconf, tmp_oconf)
			oconf = unique(oconf)
		end
		ProgressMeter.next!(p; showvalues = [
			(:Tstart, Trange[1]),
			(:T, round(T, digits=4)),
			(:Tend, round(Trange[end], digits=4)),
			(:Fmin, Fmin)
		])
	end
	return oconf,F
end

"""
	count_mismatches(g::Graph)

Count the number of topological mismatches in `g`.
  Equivalent to `compute_energy(conf, g)` with `conf = ones(Bool)`.
"""
function count_mismatches(g::Graph)
	conf = ones(Bool, length(g.leaves))
	return compute_energy(conf, g)
end

"""
	count_mismatches(t::Vararg{Tree})
"""
function count_mismatches(trees::Vararg{Tree})
	treelist = [copy(convert(Tree{TreeTools.EmptyData}, t)) for t in trees]
	mcc = naive_mccs(treelist)
	mcc_names = TreeKnit.name_mcc_clades!(treelist, mcc)
	for (i,t) in enumerate(treelist)
		TreeKnit.reduce_to_mcc!(t, mcc)
	end
	g = trees2graph(treelist)
	conf = ones(Bool, length(g.leaves))

	return compute_energy(conf, g)
end


