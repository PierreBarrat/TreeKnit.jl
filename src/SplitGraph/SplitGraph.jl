module SplitGraph

using TreeKnit
using TreeTools
using ProgressMeter
using LoggingExtras

include("objects.jl")
include("tools.jl")
include("energy.jl")
include("likelihood.jl")

"""
	opttrees!(t... ; kwargs...)

Return a list of MCCs for input trees.
Output:
1.
"""
function opttrees(
	t...;
	γ = 2,
	seq_lengths = 1000 * ones(Int64, length(t)),
	Trange = reverse(0.001:0.01:1.),
	M = 10,
	likelihood_sort = true,
	resolve = true,
	sa_rep = 1,
)
	opttrees!(
		γ, Trange, M, seq_lengths, [copy(x) for x in t]...;
		likelihood_sort, resolve, sa_rep
	)
end

function opttrees!(
	γ, Trange, M, seq_lengths, t::Vararg{Tree}; 
	likelihood_sort=true, resolve=true, sa_rep=1
)

	treelist = t
	mcc = naive_mccs(treelist...)
	if length(mcc) == 1
		return mcc, 0, 0., Int64[], Float64[], Union{Missing,Float64}[]
	end
	mcc_names = TreeKnit.name_mcc_clades!(treelist, mcc)
	for (i,t) in enumerate(treelist)
		TreeKnit.reduce_to_mcc!(t, mcc)
	end
	g = trees2graph(treelist)

	# SA - Optimization
	oconfs, F, nfound = sa_opt(g; Trange, γ, M, rep=sa_rep, resolve)
	# Computing likelihoods
	if length(oconfs) != 1
		@logmsg LogLevel(-1) "Sorting $(length(oconfs)) topologically equivalent configurations."
		@logmsg LogLevel(-2) "Configurations\n $oconfs"
		@logmsg LogLevel(-2) g.labels
		oconf, L = sortconf(oconfs, treelist, g, seq_lengths, mcc_names, likelihood_sort)
	else
		oconf = oconfs[1]
		L = Union{Missing,Float64}[]
	end
	@logmsg LogLevel(-2) "Final configuration for this iteration: $oconf."
	@logmsg LogLevel(-2) "MCCs removed: $([mcc_names[x] for x in g.labels[.!oconf]])"
	return (
		[mcc_names[x] for x in g.labels[.!oconf]],
		compute_energy(oconf,g),
		compute_F(oconf, g, γ),
		L
	)
end


function sortconf(oconfs, trees, g, seq_lengths, mcc_names, likelihood_sort, E_sort=false)
	if E_sort # Only considering configurations of lowest energies
		E = [compute_energy(conf,g) for conf in oconfs]
		Emin = minimum(E)
		oconfs_ = oconfs[findall(x->x==Emin, E)]
		@logmsg LogLevel(-2) "Removing ", length(oconfs) - length(oconfs_), " configurations using energy."
	else # Removing configurations where nothing is removed
		oconfs_ = oconfs[findall(c->sum(c)<length(c), oconfs)]
	end

	# Sorting the remaining configurations using Poisson likelihood ratios
	if length(oconfs_) == 1
		return oconfs_[1], Union{Missing,Float64}[]
	elseif !likelihood_sort
		@logmsg LogLevel(-2) "No likelihood sort: picking a random configuration among remaining ones"
		return rand(oconfs_), Union{Missing,Float64}[]
	else
		@logmsg LogLevel(-2) "Comparing $(length(oconfs_)) configurations using likelihood"
		L = Union{Missing,Float64}[]
		for conf in oconfs_
			push!(L, conf_likelihood(conf, g, seq_lengths, trees, mode=:time))
		end
		@logmsg LogLevel(-2) configurations=[[mcc_names[x] for x in g.labels[.!conf]] for conf in oconfs_]
		@logmsg LogLevel(-2) "Likelihoods: $L"
		Lmax = maximum(L)
		ismissing(Lmax) && @warn "Maximum likelihood is `missing`: this may be due to missing branch lengths"
		oconfs_ = oconfs_[findall(isequal(Lmax), L)]
		if length(oconfs_) != 1 # Final sort by energy if more than one most likely conf
			E = [compute_energy(conf,g) for conf in oconfs_]
			Emin = minimum(E)
			oconfs_ = oconfs_[findall(isequal(Emin), E)]
		end
		return rand(oconfs_), L
	end
end

end #module
