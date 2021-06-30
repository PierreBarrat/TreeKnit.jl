"""
Core function for computing splits in MCCs. Include root of MCC if it's not root of tree.
"""
function _splits_in_mcc(
	m::Array{<:AbstractString},
	t::Tree,
	leaves::Array{<:AbstractString},
	leafmap::Dict
)
    # Mask corresponding to leaves in the MCC
    mask = zeros(Bool, length(leaves))
    for l in m
        mask[leafmap[l]] = true
    end
    #
    r = lca(t.lleaves[x] for x in m)
    S = SplitList(r, leaves, mask, leafmap)
    TreeTools.clean!(S, clean_root=false)
    return S
end
"""
    splits_in_mcc(m::Array{<:AbstractString}, t::Tree)

Compute splits in MCC `m` in tree `t`. Return a `SplitList` object, with `mask` corresponding to leaves in `m`.
"""
function splits_in_mcc(m::Array{<:AbstractString}, t::Tree)
    leaves = sort(collect(keys(t.lleaves)))
    leafmap = Dict(leaf=>i for (i,leaf) in enumerate(leaves))
    return _splits_in_mcc(m, t, leaves, leafmap)
end
splits_in_mcc(m::Array{<:AbstractString}, trees::Vararg{Tree}) = Tuple(splits_in_mcc(m,t) for t in trees)

"""
    splits_in_mccs(MCCs, t::Vararg{Tree})

Compute splits in `MCCs` in tree `t`. Return an array of `SplitList` objects, with `mask` corresponding to leaves in each MCC.
When called with multiple trees, return a tuple of arrays of `SplitList` objects, each of which corresponds to one of the trees (order conserved).
"""
function splits_in_mccs(MCCs, t::Tree)
    leaves = sort(collect(keys(t.lleaves)))
    leafmap = Dict(leaf=>i for (i,leaf) in enumerate(leaves))
    return [_splits_in_mcc(m, t, leaves, leafmap) for m in MCCs]
end
splits_in_mccs(MCCs, trees::Vararg{Tree}) = Tuple(splits_in_mccs(MCCs,t) for t in trees)


"""
    new_splits(MCCs, t1::Tree, t2::Tree)

Given MCCs and two trees, what are the new splits introduced in either of the trees.
Return an array `S` of `SplitList` objects, with `S[i][m]` corresponding to new splits in tree `i` from MCC `m`.

*Note*: Resulting splits are not always unique once constrained to leaves of MCCs!
"""
function new_splits(MCCs, t1::Tree, t2::Tree)
    splits1 = SplitList(t1)
    splits2 = SplitList(t2)
    # Splits corresponding to each mcc in both trees
    MCC_splits = splits_in_mccs(MCCs, t1, t2)
    MCC_splits[1] = TreeTools.map_splits_to_tree(MCC_splits[1], t2)
    MCC_splits[2] = TreeTools.map_splits_to_tree(MCC_splits[2], t1)
    # Are those new or not?
    _new_splits!(MCC_splits[2], splits1)
    _new_splits!(MCC_splits[1], splits2)
    #
    return MCC_splits[end:-1:1]
end

"""
    new_splits(trees::Dict, MCCs::Dict)

Find new splits in each tree of `trees` according to other trees and to `MCCs`. `MCCs` should be a dictionary indexed by pairs of keys of `trees`.
Return a dictionary of the form `ns[s]` with `s` a key of `trees`.
`ns[s]` is an array of `SplitList` objects, mapped onto `trees[s]`, with length `length(trees)-1`.
"""
function new_splits(trees::Dict, MCCs::Dict)
    ns = Dict{Any, Array{SplitList,1}}()
    for (sref,tref) in trees
        ns[sref] = Array{SplitList,1}(undef, length(trees)-1)
        i = 1
        for (s,t) in trees
            if s != sref
                ns[sref][i] = RecombTools.new_splits(tref, MCCs[sref,s], t)
                i += 1
            end
        end
    end
    return ns
end

"""
    new_splits(tref::Tree, MCCs, t::Tree)

What are the splits introduced in `tref` from `t` using consistent clades `MCCs`?
Return a single `SplitList` object, with unique splits mapped onto the leaves of `tref`.
"""
function new_splits(tref::Tree, MCCs, t::Tree)
    # Splits in `tref`
    S_ref = SplitList(tref)
    # Splits corresponding to each mcc in tree `t`
    MCC_splits = TreeTools.map_splits_to_tree(splits_in_mccs(MCCs, t), tref)
    # Take the new ones only
    _new_splits!(MCC_splits, S_ref)
    # Map them onto leaves of `tref`
    return unique(TreeTools.map_splits_to_tree(MCC_splits, tref), usemask=false)
end

function _new_splits!(MCC_splits::SplitList, tree_splits::SplitList)
    idx = Int64[]
    for (i,s) in enumerate(MCC_splits)
        in(s, tree_splits, MCC_splits.mask) && push!(idx,i)
    end
    deleteat!(MCC_splits.splits, idx)
end


#=
To move to `artificialdata.jl`?
=#
"""
Indices of true splits in `S` w.r. to `Sref`.
"""
function true_splits(S::SplitList, Sref::SplitList, mask=S.mask)
    idx = Int64[]
    for (i,s) in enumerate(S)
        in(s, Sref, mask) && push!(idx, i)
    end
    return idx
end
true_splits(S, tref::Tree, mask=S.mask) = true_splits(S, SplitList(tref), mask)
