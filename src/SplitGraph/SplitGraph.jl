module SplitGraph

using RecombTools
using TreeTools

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


opttrees(t...; kwargs...) = opttrees!([copy(x, TreeTools.EmptyData) for x in t]...; kwargs...)
function opttrees(γ, Trange, M, seq_lengths, t...)
	opttrees!(γ, Trange, M, seq_lengths, [copy(x, TreeTools.EmptyData) for x in t]...)
end
"""
	opttrees!(t... ; kwargs...)

Return a list of MCCs for input trees.
Output:
1.
"""
function opttrees!(t...;
	γ=1.05,
	seq_lengths=1000 * ones(Int64, length(t)),
	Trange=reverse(0.01:0.05:1.1),
	M = 1000,
	likelihood_sort=true,
	resolve=true,
	sa_rep = 1
)
	opttrees!(
		γ, Trange, M, seq_lengths, t...;
		likelihood_sort, resolve, sa_rep,
	)
end
function opttrees!(
	γ, Trange, M, seq_lengths, t::Vararg{Tree};
	likelihood_sort=true, resolve=true, sa_rep=1,
)
	treelist = convert(Vector{Any}, collect(t))
	mcc = naive_mccs(treelist)
	if length(mcc) == 1
		return mcc, 0, 0., Int64[], Float64[], Union{Missing,Float64}[]
	end
	mcc_names = RecombTools.name_mcc_clades!(treelist, mcc)
	for (i,t) in enumerate(treelist)
		treelist[i] = RecombTools.reduce_to_mcc(t, mcc)
	end
	g = trees2graph(treelist)

	# SA - Optimization
	oconfs, E, F, nfound = sa_opt(g, γ=γ, Trange=Trange, M=M, rep=sa_rep, resolve=resolve)
	# Computing likelihoods
	if length(oconfs) != 1
		v() && println("Sorting $(length(oconfs)) topologically equivalent configurations.")
		vv() && println("Configurations\n", oconfs)
		vv() && println(g.labels)
		oconf, L = sortconf(oconfs, treelist, g, seq_lengths, mcc_names, likelihood_sort, false)
	else
		oconf = oconfs[1]
		L = Union{Missing,Float64}[]
	end
	vv() && println("Final configuration for this iteration: $oconf.")
	vv() && println("MCCs removed: $([mcc_names[x] for x in g.labels[.!oconf]])")
	return [mcc_names[x] for x in g.labels[.!oconf]], compute_energy(oconf,g), compute_F(oconf, g, γ), E, F, L
end


function sortconf(oconfs, trees, g::Graph, seq_lengths, mcc_names, likelihood_sort, E_sort)
	if E_sort # Only considering configurations of lowest energies
		E = [compute_energy(conf,g) for conf in oconfs]
		Emin = minimum(E)
		oconfs_ = oconfs[findall(x->x==Emin, E)]
		v() && println("Removing ", length(oconfs) - length(oconfs_), " configurations using energy.")
	else # Removing configurations where nothing is removed
		oconfs_ = oconfs[findall(c->sum(c)<length(c), oconfs)]
	end

	# Sorting the remaining configurations using Poisson likelihood ratios
	if length(oconfs_) == 1
		return oconfs_[1], Union{Missing,Float64}[]
	elseif !likelihood_sort
		v() && println("Picking a random configuration among remaining ones")
		return rand(oconfs_), Union{Missing,Float64}[]
	else
		v() && println("Comparing $(length(oconfs_)) configurations using likelihood")
		L = Union{Missing,Float64}[]
		for conf in oconfs_
			vv() && println("## Looking at configuration $conf with energy $(compute_energy(conf,g))")
			push!(L, conf_likelihood(conf, g, seq_lengths, trees, mode=:time, v=vv()))
			vv() && println()
		end
		vv() && println("Confs: ", [[mcc_names[x] for x in g.labels[.!conf]] for conf in oconfs_])
		v() && println("Likelihoods: ", L)
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
