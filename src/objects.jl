####################################
### TreeKnit run options objects ###
####################################

"""
	struct OptArgs

Storing parameters for `SplitGraph.runopt` function.

### General
- `γ::Real = 2`
- `itmax::Int64 = 15`: maximal number of iterations of naive MCCs / SA cycles.
- `likelihood_sort::Bool = true`: sort equivalent configurations using likelihood test
  based on branch length.
- `resolve::Bool = true`: try to resolve trees while finding MCCs.
- `strict::Bool = true`: only resolve unambiguous splits
- `seq_lengths`: lengths of sequences that trees were built from.
  Used in likelihood calculations.
  This is initialized from other input arguments, and defaults to sequences of length one.
### Pipeline options - necessary for K>2 trees
- `pre_resolve::Bool = true`: pre-resolve all trees using each other prior to MCC inference
- `rounds::Int=1`: #rounds of MCC inference and tree resolution on all tree pairs
- `final_no_resolve::Bool = false`: do not resolve in the final round of pair-wise MCC inference
- `parallel::Bool = false`: parallelize MCC inference as much as possible. 
### Simulated annealing
- `nMCMC::Int = 25`: The *total* number of MCMC steps (swaps) for a tree of `n` leaves is `nMCMC*n`. 
  The number of MCMC steps at one temperature is `nMCMC * n / nT`.
- `cooling_schedule = :geometric`: type of cooling schedule `(:geometric, :linear, :acos)`
- `Tmin::Float64 = 0.05`: minimal temperature of SA.
- `Tmax::Float64 = 0.8`: maximal temperature of SA.
- `nT::Int = 3000`: number of steps in the cooling schedule
"""
@with_kw mutable struct OptArgs
	γ::Real = 2
	itmax::Int64 = 15
	likelihood_sort::Bool = true
	resolve::Bool = true
	strict::Bool = true 
	seq_lengths::Vector{Int} = [1, 1]
	# Pipeline options - needed for K>2 trees
	pre_resolve::Bool = true
	rounds::Int=1
	final_no_resolve::Bool = false
	parallel::Bool = false
	# For the annealing
	nMCMC::Int = 50
	sa_rep::Int64 = 1
	Tmin::Float64 = 0.05; @assert Tmin > 0
	Tmax::Float64 = 1; @assert Tmax > Tmin
	nT::Int = 100
	cooling_schedule = :geometric
	Trange = get_cooling_schedule(Tmin, Tmax, nT, type=cooling_schedule)
end


"""
	function OptArgs(K::Int; method=:None, kwargs...)

OptArgs constructor will default to :better_trees method if K>2, :better_MCCs if K==2.
"""
function OptArgs(K::Int; method=:None, kwargs...)
	kwargs_dict = Dict(kwargs)
	oa = OptArgs() #default values for K==2
    oa.seq_lengths = get(kwargs_dict, :seq_lenghts, repeat([1], K))
	oa.pre_resolve = get(kwargs_dict, :pre_resolve, true)
	oa.strict = get(kwargs_dict, :strict, true)

	if (method == :better_trees) || ((method == :None) && K>2)
        oa.resolve = get(kwargs_dict, :resolve, false)
        oa.final_no_resolve = get(kwargs_dict, :final_no_resolve, true)
		oa.rounds = get(kwargs_dict, :rounds, 1)
	elseif (method == :better_MCCs) || ((method == :None) && K==2)
        oa.resolve = get(kwargs_dict, :resolve, true)
		if K>2
            oa.final_no_resolve = get(kwargs_dict, :final_no_resolve, true)
			oa.rounds = get(kwargs_dict, :rounds, 2)
		else
            oa.final_no_resolve = get(kwargs_dict, :final_no_resolve, false)
            oa.rounds = get(kwargs_dict, :rounds, 1)
		end
	else
		error("`method` should be one of `:BetterTrees`, `:BetterMCCs` or `:None` - got $method.")
	end

	##clean up
	if (oa.final_no_resolve==true) && (oa.rounds==1) && (oa.resolve==true)
		oa.resolve=false
		oa.final_no_resolve=false
	end

	return oa
end

function get_cooling_schedule(Tmin, Tmax, nT; type=:geometric)
	if type == :geometric
		return get_geometric_cooling_schedule(Tmin, Tmax, nT)
	elseif type == :linear
		return get_linear_cooling_schedule(Tmin, Tmax, nT)
	elseif type == :acos
		return get_acos_cooling_schedule(Tmin, Tmax, nT)
	else
		error("Unknown `cooling_schedule` field: $(type). See `?OptArgs` for allowed values.")
	end
end

function get_geometric_cooling_schedule(Tmin, Tmax, nT)
	α = exp((log(Tmin) - log(Tmax)) / nT)
	n = ceil(Int, (log(Tmin) - log(Tmax)) / log(α))
	return [α^i * Tmax for i in 0:n]
end

function get_acos_cooling_schedule(Tmin, Tmax, nT)
	f(x, K=1.5) = if x < 0.5
		0.5 + 2^(K-1)*abs(acos(2*x-1)/3.14-0.5)^K
	elseif x == 0.5
		0.5
	else
		0.5 - 2^(K-1)*abs(acos(2*x-1)/3.14-0.5)^K
	end
	# f(0) = 1.0 and f(1) = 0.
	return [(Tmax - Tmin)*f(t) + Tmin for t in range(0, stop = 1, length = nT)]
end

get_linear_cooling_schedule(Tmin, Tmax, nT) = return collect(reverse(range(Tmin, stop = Tmax, length = nT)))


####################
### MCC objects ###
####################
"""
	struct MCC_set

structure to store and access computed MCCs for tree pairs

- `no_trees`: number of trees
- `order_trees`: vector of tree labels.
- `mccs`: Dictionary of calculated mccs, key is set of labels of trees in each tree pair

Add and retrieve mccs with `get!` and `add!` and Tuple or Vararg of tree labels or position of tree labels in `order_trees`
"""
struct MCC_set
    no_trees :: Int
    order_trees :: Vector{String}
    mccs :: Dict{Set{String}, Vector{Vector{String}}}
end

function MCC_set(no_trees, order_trees; mccs = Dict())
	return MCC_set(no_trees, order_trees, mccs)
end
 
function MCC_set(no_trees, order_trees, MCC_list::Vector{Vector{Vector{String}}})
    @assert length(MCC_list) >= (no_trees*(no_trees-1)/2)
    if !allunique(order_trees)
    	@warn "Tree labels are not unique: this will cause issues when labeling MCC pairs." order_trees
    end

    M = MCC_set(no_trees, order_trees)
    # iters = Combinatorics.combinations(1:no_trees, 2)
    pair_iter = ordered_pairs(no_trees)
    for (i, comb) in enumerate(pair_iter)
        add!(M, MCC_list[i], comb...)
    end
    return M
end

function Base.get(M::MCC_set, pos::Vararg{String})
    if haskey(M.mccs, Set(pos))
        return M.mccs[Set(pos)]
    else
        return nothing
    end
end

Base.get(M::MCC_set, pos::Tuple{Vararg{String}}) = get(M, pos...)
Base.get(M::MCC_set, pos::Vararg{Int}) = get(M, [M.order_trees[i] for i in pos]...)
Base.get(M::MCC_set, pos::Tuple{Vararg{Int}}) = get(M, [M.order_trees[i] for i in pos]...)

Base.getindex(M::MCC_set, pos...) = get(M, pos...) # for M[i,j] syntax

function add!(M::MCC_set, value::Vector{Vector{String}}, pos::Vararg{String})
    M.mccs[Set(pos)] = value
end

function add!(M::MCC_set, value::Vector{Vector{String}}, pos::Tuple{Vararg{String}})
    M.mccs[Set(pos)] = value
end

function add!(M::MCC_set, value::Vector{Vector{String}}, pos::Tuple{Vararg{Int}})
    add!(M, value, [M.order_trees[i] for i in pos]...)
end

function add!(M::MCC_set, value::Vector{Vector{String}}, pos::Vararg{Int})
    add!(M, value, [M.order_trees[i] for i in pos]...)
end

"""
	iter_pairs(M::MCC_set)

Iterate over all tree pairs and their calculated MCCs in order of calculation.
"""
function iter_pairs(M::MCC_set)
    pair_iter = ordered_pairs(M.no_trees)
    return [M.order_trees[p] for p in pair_iter], [get(M, p...) for p in pair_iter]
end

"""
	iter_pairs(M::MCC_set, tree_label::String)

Iterate over all tree pairs with tree `tree_label` and their calculated MCCs in order of calculation.
"""
function iter_shared(M::MCC_set, tree_label::String)
    return [get(M, (tree, t)) for t in M.tree_order if t!=tree]
end

function Base.print(M::MCC_set)
    for (key, value) in M.mccs
        println("")
        println(key)
        println(value)
    end
end

function Base.copy(M::MCC_set)
    copy_M = MCC_set(M.no_trees, deepcopy(M.order_trees); mccs =deepcopy(M.mccs))
    return copy_M
end

ordered_pairs(K) = (([i,j] for j in (i+1):K) for i in 1:(K-1)) |> Iterators.flatten

