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
function opttrees(treelist::Vector{Tree{T}}, copyleaves::Union{Nothing, Dict{Vector{Int64}, Set{String}}};
	mcc_names=Dict(),
	γ=2,
	seq_lengths=1000 * ones(Int64, length(treelist)),
	Trange=reverse(0.001:0.01:1.),
	M = 10,
	likelihood_sort=true,
	resolve=true,
	sa_rep = 1,
	verbose=true
) where T 
	opttrees!(
		γ, Trange, M, seq_lengths, [copy(t) for t in treelist], deepcopy(copyleaves);
		mcc_names, likelihood_sort, resolve, sa_rep, verbose
	)
end
function opttrees!(
	γ, Trange, M, seq_lengths, treelist::Vector{Tree{T}}, copyleaves; mcc_names=Dict(),
	likelihood_sort=true, resolve=true, sa_rep=1, verbose=true,
) where T 
	set_verbose(verbose)
	if isnothing(copyleaves)
		print("no copy leaves given")
		treelist, copyleaves = TreeKnit.prepare_copies!(treelist);
	end
	mcc = Dict()
	for tree_pair in Combinatorics.combinations(1:length(treelist), 2)
		mcc[tree_pair] = naive_mccs([treelist[tree_pair[1]], treelist[tree_pair[2]]], copyleaves[tree_pair])
	end
	mcc_names = TreeKnit.name_mcc_clades!(treelist, copyleaves, mcc; mcc_names)
	for t in treelist
		t = TreeKnit.remove_zero_copies!(t)
	end
	if [length(mcc[pair]) for pair in Combinatorics.combinations(1:length(treelist), 2)] == ones(Int(length(treelist)*(length(treelist)-1)/2))
		return [mcc_names[[1,2]][x] for x in mcc[[1,2]]], 0, 0., Int64[], Float64[], Union{Missing,Float64}[]
	end
	print("after naive")
	print(copyleaves)
	print(treelist[1])
	print(treelist[2])
	g = trees2graph(treelist)
	# SA - Optimization
	oconfs, F, nfound = sa_opt(g, γ=γ, Trange=Trange, M=M, rep=sa_rep, resolve=resolve)
	# Computing likelihoods
	if length(oconfs) != 1
		v() && @info "Sorting $(length(oconfs)) topologically equivalent configurations."
		vv() && @info "Configurations\n $oconfs"
		vv() && @info g.labels
		oconf, L = sortconf(oconfs, treelist, g, seq_lengths, mcc_names[[1,2]], likelihood_sort, false)
	else
		oconf = oconfs[1]
		L = Union{Missing,Float64}[]
	end
	vv() && @info "Final configuration for this iteration: $oconf."
	vv() && @info "MCCs removed: $([mcc_names[[1,2]][x] for x in g.labels[.!oconf]])"
	return [mcc_names[[1,2]][x] for x in g.labels[.!oconf]], compute_energy(oconf,g), compute_F(oconf, g, γ), L
end

function opttrees(trees::Vararg{Tree};
	γ=2,
	seq_lengths=1000 * ones(Int64, length(trees)),
	Trange=reverse(0.001:0.01:1.),
	M = 10,
	likelihood_sort=true,
	resolve=true,
	sa_rep = 1,
	verbose=false
) where T 
return opttrees([copy(t) for t in trees], nothing; γ, seq_lengths, Trange, M, likelihood_sort, resolve, sa_rep, verbose)
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



end # module
