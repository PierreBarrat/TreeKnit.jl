export MCC_set 

struct MCC_set
    no_trees :: Int
    order_trees :: Vector{String}
    mccs :: Dict{Set{String}, Vector{Vector{String}}}
end

MCC_set(no_trees::Int, order_trees:: Vector{String}; mccs =Dict{Set{String}, Vector{Vector{String}}}()) = MCC_set(no_trees, order_trees, mccs)
 
function Base.get(M::MCC_set, pos::Vararg{String})
    if haskey(M.mccs, Set(pos))
        return M.mccs[Set(pos)]
    else
        return nothing
    end
end

function Base.get(M::MCC_set, pos::Tuple{Vararg{String}})
    return get(M, pos...)
end

function Base.get(M::MCC_set, pos::Vararg{Int})
    return get(M, [M.order_trees[i] for i in pos]...)
end

function Base.get(M::MCC_set, pos::Tuple{Vararg{Int}})
    return get(M, [M.order_trees[i] for i in pos]...)
end

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

function iter_pairs(M::MCC_set)
    pairs = Combinatorics.combinations(1:M.no_trees, 2)
    return [M.order_trees[p] for p in pairs], [get(M, p...) for p in pairs]
end

function iter_shared(M::MCC_set, tree::String)
    return [get(M, (tree, t)) for t in M.tree_order if t!=tree]
end

function convert_MCC_list_to_set(no_trees::Int, order_trees::Vector{String}, MCC_list::Vector{Vector{Vector{String}}})
    @assert length(MCC_list) >= (no_trees*(no_trees-1)/2)
    M = MCC_set(no_trees, order_trees)
    iters = Combinatorics.combinations(1:no_trees, 2)
    for (i, comb) in enumerate(iters)
        TreeKnit.add!(M, MCC_list[i], comb...)
    end
    return M
end

function Base.print(M::MCC_set)
    for (key, value) in M.mccs
        println("")
        println(key)
        println(value)
    end
end