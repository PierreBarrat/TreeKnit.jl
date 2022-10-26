module SplitGraph

using TreeKnit
using TreeTools
using ProgressMeter

include("objects.jl")
include("tools.jl")
include("energy.jl")
include("likelihood.jl")

let verbose::Bool = false, vverbose::Bool = false
	global v() = verbose
	global vv() = vverbose
	global set_verbose(v) = (verbose = v)
	global set_vverbose(v) = (vverbose = v)
end


# function opttrees(γ, Trange, M, seq_lengths, t...)
# 	opttrees!(γ, Trange, M, seq_lengths, [copy(x, TreeTools.EmptyData) for x in t]...)
# end

# function opttrees(t...; kwargs...)
# 	opttrees!([copy(x, TreeTools.EmptyData) for x in t]...; kwargs...)
# end
"""
	opttrees!(t... ; kwargs...)

Return a list of MCCs for input trees.
Output:
1.
"""
function opttrees(t...;
	γ=2,
	seq_lengths=1000 * ones(Int64, length(t)),
	Trange=reverse(0.001:0.01:1.),
	M = 10,
	likelihood_sort=true,
	resolve=true,
	sa_rep = 1,
	consistent =false,
	constraint_cost= γ,
	verbose=false,
	shared_maps=nothing
)
	opttrees!(
		γ, Trange, M, seq_lengths, [copy(convert(Tree{TreeTools.MiscData}, x)) for x in t]...;
		likelihood_sort=likelihood_sort, resolve=resolve, sa_rep=sa_rep, consistent=consistent, 
		constraint_cost=constraint_cost, verbose=verbose, shared_maps=shared_maps
	)
end

function opttrees!(
	γ, Trange, M, seq_lengths, t::Vararg{Tree}; 
	likelihood_sort=true, resolve=true, sa_rep=1, consistent =false, constraint_cost= γ, verbose=false, shared_maps=nothing
)
	set_verbose(verbose)

	treelist = t
	mcc = naive_mccs(treelist...)
	if length(mcc) == 1
		return mcc, 0, 0., Int64[], Float64[], Union{Missing,Float64}[]
	end
	mcc_names = TreeKnit.name_mcc_clades!(treelist, mcc; shared_maps)
	for (i,t) in enumerate(treelist)
		TreeKnit.reduce_to_mcc!(t, mcc)
	end
	g = trees2graph(treelist)
	if consistent
		mask = get_consistency_mask(g, t...; shared_maps)
	else
		mask =  []
	end

	# SA - Optimization
	oconfs, F, nfound = sa_opt(g, γ=γ, Trange=Trange, M=M, rep=sa_rep, resolve=resolve, mask=mask, constraint_cost=constraint_cost)
	# Computing likelihoods
	if length(oconfs) != 1
		v() && @info "Sorting $(length(oconfs)) topologically equivalent configurations."
		vv() && @info "Configurations\n $oconfs"
		vv() && @info g.labels
		oconf, L = sortconf(oconfs, treelist, g, seq_lengths, mcc_names, likelihood_sort, false)
	else
		oconf = oconfs[1]
		L = Union{Missing,Float64}[]
	end
	vv() && @info "Final configuration for this iteration: $oconf."
	vv() && @info "MCCs removed: $([mcc_names[x] for x in g.labels[.!oconf]])"
	return [mcc_names[x] for x in g.labels[.!oconf]], compute_energy(oconf,g), compute_F(oconf, g, γ, mask=mask, constraint_cost=constraint_cost), L
end


function sortconf(oconfs, trees, g::Graph, seq_lengths, mcc_names, likelihood_sort, E_sort)
	if E_sort # Only considering configurations of lowest energies
		E = [compute_energy(conf,g) for conf in oconfs]
		Emin = minimum(E)
		oconfs_ = oconfs[findall(x->x==Emin, E)]
		v() && @info "Removing ", length(oconfs) - length(oconfs_), " configurations using energy."
	else # Removing configurations where nothing is removed
		oconfs_ = oconfs[findall(c->sum(c)<length(c), oconfs)]
	end

	# Sorting the remaining configurations using Poisson likelihood ratios
	if length(oconfs_) == 1
		return oconfs_[1], Union{Missing,Float64}[]
	elseif !likelihood_sort
		v() && @info "No likelihood sort: picking a random configuration among remaining ones"
		return rand(oconfs_), Union{Missing,Float64}[]
	else
		v() && @info "Comparing $(length(oconfs_)) configurations using likelihood"
		L = Union{Missing,Float64}[]
		for conf in oconfs_
			push!(L, conf_likelihood(conf, g, seq_lengths, trees, mode=:time, v=vv()))
		end
		vv() && println("Confs: ", [[mcc_names[x] for x in g.labels[.!conf]] for conf in oconfs_])
		vv() && @info "Likelihoods: $L"
		Lmax = maximum(L)
		ismissing(Lmax) && @warn "Maximum likelihood is `missing`"
		oconfs_ = oconfs_[findall(isequal(Lmax), L)]
		if length(oconfs_) != 1 # Final sort by energy if more than one most likely conf
			E = [compute_energy(conf,g) for conf in oconfs_]
			Emin = minimum(E)
			oconfs_ = oconfs_[findall(isequal(Emin), E)]
		end
		return rand(oconfs_), L
	end
end

"""
	get_consistency_mask(g, tree)

Get which branches/ terminal nodes of the current tree (SplitGraph) 
should not be removed in the subsequent MCMC due to the fact that 
they have a `shared_branch`.
"""
function get_consistency_mask(g::Graph, trees::Vararg{Tree}; shared_maps=nothing)
	if isnothing(shared_maps)
		return []
	end

	mask_names = String[]
	for (i, tree) in enumerate(trees)
		for leaf in POTleaves(tree)
			if shared_maps[i][leaf.label]
				push!(mask_names, leaf.label)
			end
		end
	end
	mask_names = sort!(unique!(mask_names))

	mask = [g.labels_to_int[m] for m in mask_names]

	return mask
end

end # module
