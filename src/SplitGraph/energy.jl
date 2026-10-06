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
					elseif !get_resolve() && !are_equal(a1, a2, conf)
						E += 1
					end
				end
			end
		end
	end
	return E
end

"""
	are_equal(a1::SplitNode, a2::SplitNode, conf)

Are the splits of `a1` and `a2` equal for leaf configuration `conf`?
"""
function are_equal(a1::SplitNode, a2::SplitNode, conf)
	return is_contained(a1, a2, conf) && is_contained(a2, a1, conf)
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
	if are_equal(a1, a2, conf)
		return true
	elseif is_contained(a1, a2, conf) && _are_equal_with_resolution(g, a1, a2, conf)
		return true
	elseif is_contained(a2, a1, conf) && _are_equal_with_resolution(g, a2, a1, conf)
		return true
	end
	return false
end

"""
	_are_equal_with_resolution(g, a1, a2, conf)

For every internal child of `a2`, check that its split is either disjoint from the split of
`a1` or included in it.
"""
function _are_equal_with_resolution(g::SplitGraph.Graph, a1::SplitNode, a2::SplitNode, conf)
	for c in a2.child
		if c isa SplitNode && !is_contained(c, a1, conf) && !are_disjoint(c, a1, conf)
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
			n >= 2 && return false
		end
	end
	return true
end

"""
	is_contained(a1::SplitNode, a2::SplitNode, conf)

Is the split of `a1` included in the split of `a2` for leaf configuration `conf`?
"""
function is_contained(a1::SplitNode, a2::SplitNode, conf)
	@inbounds for i in a1.conf
		if conf[i] && !a2.mask[i]
			return false
		end
	end
	return true
end

"""
	are_disjoint(a1::SplitNode, a2::SplitNode, conf)

Are the splits of `a1` and `a2` disjoint for leaf configuration `conf`?
"""
function are_disjoint(a1::SplitNode, a2::SplitNode, conf)
	small, big = length(a1.conf) > length(a2.conf) ? (a2, a1) : (a1, a2)
	@inbounds for i in small.conf
		if conf[i] && big.mask[i]
			return false
		end
	end
	return true
end


"""
	SplitLayout(g::Graph)

Integer layout of the internal nodes of `g` for `EnergyState`, numbered as in `g.internals`.
- `parent[a]`: parent of node `a`, `0` for a root
- `splitchild[a]`, `leafchild[a]`: internal children and leaf children of `a`
- `clade[:, a]`: the leaves below `a` as a bit set in 64-bit words
- `leafanc[k, i]`: parent of leaf `i` in tree `k`
"""
struct SplitLayout
	K::Int
	trivial::Bool
	parent::Vector{Int}
	splitchild::Vector{Vector{Int}}
	leafchild::Vector{Vector{Int}}
	clade::Matrix{UInt64}
	leafanc::Matrix{Int}
end

function SplitLayout(g::Graph)
	index = IdDict{SplitNode,Int}(a => n for (n, a) in enumerate(g.internals))
	nnodes = length(g.internals)
	nleaves = length(g.leaves)
	parent = zeros(Int, nnodes)
	splitchild = [Int[] for _ in 1:nnodes]
	leafchild = [Int[] for _ in 1:nnodes]
	clade = zeros(UInt64, cld(nleaves, 64), nnodes)
	for (a, node) in enumerate(g.internals)
		node.isroot || (parent[a] = index[node.anc])
		for c in node.child
			if c isa SplitNode
				push!(splitchild[a], index[c])
			else
				push!(leafchild[a], c.index)
			end
		end
		for i in node.conf
			clade[(i - 1) >> 6 + 1, a] |= one(UInt64) << ((i - 1) & 63)
		end
	end
	leafanc = [index[g.leaves[i].anc[k]] for k in 1:g.K, i in 1:nleaves]
	trivial = nleaves + nnodes == 1
	return SplitLayout(g.K, trivial, parent, splitchild, leafchild, clade, leafanc)
end

"""
	EnergyState(L::SplitLayout, conf, resolve; terms=true)

Energy of the leaf configuration `conf`, kept equal to `compute_energy` under single-leaf
flips (`flip!`, `undo!`). With `terms=false`, only the numbers of kept leaves below the
nodes are computed, for the comparisons of `conf_likelihood`.

The energy is a sum of one term per kept leaf, which compares the first non-trivial ancestors
of the leaf in each tree. Flipping leaf `j` changes the term of `j` and of the kept leaves
whose first non-trivial ancestor in some tree is an ancestor of `j`. Each of these is the
only kept leaf below a child of an ancestor of `j`, before or after the flip. The number of
kept leaves below every node is maintained along the root paths of `j`, so these leaves and
the first non-trivial ancestors are found without scanning splits.
"""
mutable struct EnergyState
	L::SplitLayout
	resolve::Bool
	conf::Vector{Bool}
	words::Vector{UInt64}
	count::Vector{Int}
	term::Vector{Int}
	E::Int
	cand::Vector{Int}
	stamp::Vector{Int}
	step::Int
	undo::Vector{Tuple{Int,Int}}
	flipped::Int
	anc::Vector{Int}
end

function EnergyState(L::SplitLayout, conf, resolve; terms=true)
	n = length(conf)
	words = zeros(UInt64, size(L.clade, 1))
	for i in 1:n
		conf[i] && (words[(i - 1) >> 6 + 1] |= one(UInt64) << ((i - 1) & 63))
	end
	count = [sum(w -> count_ones(L.clade[w, a] & words[w]), eachindex(words); init = 0) for a in axes(L.clade, 2)]
	st = EnergyState(
		L, resolve, copy(conf), words, count, zeros(Int, n), 0, Int[], zeros(Int, n), 0, Tuple{Int,Int}[], 0,
		zeros(Int, L.K),
	)
	if terms
		for i in 1:n
			st.term[i] = leaf_term(st, i)
		end
		st.E = sum(st.term; init = 0)
	end
	return st
end

"""
	flip!(st::EnergyState, j)

Flip leaf `j` and return the new energy. `undo!` reverts the last flip.
"""
function flip!(st::EnergyState, j)
	empty!(st.cand)
	st.step += 1
	add_candidate!(st, j)
	affected_leaves!(st, j)
	toggle!(st, j)
	affected_leaves!(st, j)
	empty!(st.undo)
	@inbounds for i in st.cand
		t = leaf_term(st, i)
		push!(st.undo, (i, st.term[i]))
		st.E += t - st.term[i]
		st.term[i] = t
	end
	st.flipped = j
	return st.E
end

"""
	undo!(st::EnergyState)

Revert the last `flip!`.
"""
function undo!(st::EnergyState)
	@inbounds for (i, t) in st.undo
		st.E += t - st.term[i]
		st.term[i] = t
	end
	toggle!(st, st.flipped)
	return nothing
end

function toggle!(st::EnergyState, j)
	L = st.L
	st.conf[j] = !st.conf[j]
	st.words[(j - 1) >> 6 + 1] ⊻= one(UInt64) << ((j - 1) & 63)
	d = st.conf[j] ? 1 : -1
	@inbounds for k in 1:L.K
		a = L.leafanc[k, j]
		while a != 0
			st.count[a] += d
			a = L.parent[a]
		end
	end
	return nothing
end

# Add to `st.cand` the kept leaves that are the only kept leaf below a child of an ancestor of `j`.
function affected_leaves!(st::EnergyState, j)
	L = st.L
	@inbounds for k in 1:L.K
		a = L.leafanc[k, j]
		while a != 0
			for l in L.leafchild[a]
				st.conf[l] && add_candidate!(st, l)
			end
			for c in L.splitchild[a]
				st.count[c] == 1 && add_candidate!(st, single_kept_leaf(st, c))
			end
			a = L.parent[a]
		end
	end
	return nothing
end

# Add leaf `l` to `st.cand` once per flip.
@inline function add_candidate!(st::EnergyState, l)
	@inbounds if st.stamp[l] != st.step
		st.stamp[l] = st.step
		push!(st.cand, l)
	end
	return nothing
end

# The only kept leaf below node `a`, which has exactly one.
function single_kept_leaf(st::EnergyState, a)
	L = st.L
	@inbounds while true
		for l in L.leafchild[a]
			st.conf[l] && return l
		end
		for c in L.splitchild[a]
			if st.count[c] == 1
				a = c
				break
			end
		end
	end
end

# Contribution of leaf `i` to the energy: the number of tree pairs whose first non-trivial
# ancestors of `i` are not equal (or not compatible by resolution).
function leaf_term(st::EnergyState, i)
	L = st.L
	(st.conf[i] && !L.trivial) || return 0
	@inbounds for k in 1:L.K
		a = L.leafanc[k, i]
		while st.count[a] < 2 && L.parent[a] != 0
			a = L.parent[a]
		end
		st.anc[k] = a
	end
	E = 0
	@inbounds for k1 in 1:L.K, k2 in (k1 + 1):L.K
		compatible(st, st.anc[k1], st.anc[k2]) || (E += 1)
	end
	return E
end

# Word-level version of `are_equal` and `are_equal_with_resolution` on nodes `a1`, `a2`.
# Equal splits have equal numbers of kept leaves, and a split has no more kept leaves than a
# split that contains it, so the counts rule out most comparisons.
function compatible(st::EnergyState, a1, a2)
	c1, c2 = st.count[a1], st.count[a2]
	c1 == c2 && eq_on(st, a1, a2) && return true
	st.resolve || return false
	return (c1 <= c2 && subset_on(st, a1, a2) && refinable(st, a1, a2)) ||
		(c2 <= c1 && subset_on(st, a2, a1) && refinable(st, a2, a1))
end

# Is every internal child of `a` included in or disjoint from the split of `s`?
function refinable(st::EnergyState, s, a)
	@inbounds for c in st.L.splitchild[a]
		subset_on(st, c, s) || disjoint_on(st, c, s) || return false
	end
	return true
end

function eq_on(st::EnergyState, a, b)
	C, m = st.L.clade, st.words
	@inbounds for w in eachindex(m)
		iszero((C[w, a] ⊻ C[w, b]) & m[w]) || return false
	end
	return true
end

function subset_on(st::EnergyState, a, b)
	C, m = st.L.clade, st.words
	@inbounds for w in eachindex(m)
		iszero(C[w, a] & ~C[w, b] & m[w]) || return false
	end
	return true
end

function disjoint_on(st::EnergyState, a, b)
	C, m = st.L.clade, st.words
	@inbounds for w in eachindex(m)
		iszero(C[w, a] & C[w, b] & m[w]) || return false
	end
	return true
end

"""
"""
function compute_F(conf::Array{Bool,1}, g::Graph, γ::Real)
	E = compute_energy(conf, g)
	return E + γ*(length(conf) - sum(conf))
end

"""
	doMCMC(g::Graph, conf, M; T=1, γ=1)
	doMCMC(st::EnergyState, M; T=1, γ=1)

Run `M` Metropolis steps at temperature `T` from configuration `conf`, or from the
configuration of `st`, which the second form updates. Return the configurations of minimal
free energy, the final configuration, and the minimal free energy.
"""
function doMCMC(
	g::Graph, conf::Array{Bool,1}, M::Int64;
	T=1, γ=1
)
	return doMCMC(EnergyState(SplitLayout(g), conf, get_resolve()), M; T, γ)
end
function doMCMC(st::EnergyState, M::Int64; T=1, γ=1)
	_conf = st.conf
	nremoved = length(_conf) - sum(_conf)
	F = st.E + γ*nremoved
	Fmin = F

	oconf = [copy(_conf)]
	# `_conf` is in `oconf` and has not changed since
	listed = true
	for m in 1:M
		F, nremoved, accepted = mcmcstep!(st, F, nremoved, T, γ)
		accepted && (listed = false)
		# If new minimum is found
		if F < Fmin
			Fmin = F
			oconf[1] .= _conf
			deleteat!(oconf, 2:length(oconf))
			listed = true
		end
		# If equal minimum is found
		if F == Fmin && !listed
			in(_conf, oconf) || push!(oconf, copy(_conf))
			listed = true
		end
	end
	return oconf, _conf, Fmin
end

"""
	mcmcstep!(st::EnergyState, F, nremoved, T, γ)

Perform an mcmc step by flipping a leaf of `st` at random. `F` and `nremoved` are the free
energy and the number of removed leaves before the step.

Return `(F, nremoved, accepted)` after the step.
"""
function mcmcstep!(st::EnergyState, F, nremoved, T, γ)
	i = rand(1:length(st.conf))
	Enew = flip!(st, i)
	nremoved_new = st.conf[i] ? nremoved - 1 : nremoved + 1
	Fnew = Enew + γ*nremoved_new
	if Fnew < F || exp(-(Fnew-F)/T) > rand()
		return Fnew, nremoved_new, true
	else
		undo!(st)
		return F, nremoved, false
	end
end

"""
	sa_opt(g::Graph; Trange=1.:-0.01:0.1, γ=1.05, M=1000, rep=1, resolve=true)

Call `_sa_opt` repeatedly to find a set of optimal confs.
"""
function sa_opt(
	g::Graph;
	Trange=OptArgs().Trange, γ=2., M=10, rep=1, resolve=true
)
	set_resolve(resolve)
	#
	oconf = Any[]
	F = Int64[]
	Fmin = Inf
	nfound = 0
	for r in 1:rep
		oconf_, F_ = _sa_opt(g, γ, Trange, M)
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

function _sa_opt(g::Graph, γ, Trange, M)
	reset_chance = 0.
	conf = ones(Bool, length(g.leaves))
	layout = SplitLayout(g)
	st = EnergyState(layout, conf, get_resolve())
	oconf = [copy(conf)]
	F = Float64[Inf]
	Fmin = F[1]
	p = Progress(
		length(Trange);
		dt=1.,
		desc="Simulated annealing: ",
		enabled=true,
		barglyphs=BarGlyphs("[=> ]"),
		showspeed=true,
	)
	for T in Trange
		if rand() < reset_chance
			st = EnergyState(layout, oconf[rand(1:length(oconf))], get_resolve())
		end
		tmp_oconf, conf, fmin = doMCMC(st, M; T, γ)
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


