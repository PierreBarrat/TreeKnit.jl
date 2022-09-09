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
    TreeTools.clean!(S, clean_root=false) # if true, do not resolve based on root of MCC

    return S
end
"""
    splits_in_mcc(m::Array{<:AbstractString}, t...)

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
function new_splits(MCCs, t1::Tree, t2::Tree; strict=false)
    nS = [
    	new_splits(t1, MCCs, t2; strict=strict), # Using each tree as a reference in turn
    	new_splits(t2, MCCs, t1; strict=strict),
    ]
    return nS
end

"""
    new_splits(tref::Tree, MCCs, t::Tree)

What are the splits introduced in `tref` from `t` using consistent clades `MCCs`?
Return a single `SplitList` object, with unique splits mapped onto the leaves of `tref`.
"""
function new_splits(tref::Tree, MCCs, t::Tree; strict=false)
    # Splits in `tref`
    S_ref = SplitList(tref)
    # Splits corresponding to each mcc in tree `t`
    MCC_splits = map_splits_to_tree(splits_in_mccs(MCCs, t), tref; MCCs=MCCs, strict=strict)
    # Take the new ones only
    _new_splits!(MCC_splits, S_ref)
    # Map them onto leaves of `tref`
    return unique(map_splits_to_tree(MCC_splits, tref), usemask=false)
end

function _new_splits!(MCC_splits::SplitList, tree_splits::SplitList)
    idx = Int64[]
    for (i,s) in enumerate(MCC_splits)
        if isempty(s)
            push!(idx,i)
            continue
        end
        in(s, tree_splits, MCC_splits.mask) && push!(idx,i)
    end
    deleteat!(MCC_splits.splits, idx)
end


"""
       map_splits_to_tree(S_array::Array{<:SplitList,1}, t::Tree)

Call `map_splits_to_tree(S::SplitList, t::Tree)` for all elements of `S`.
Return a single `SplitList`.
"""
function map_splits_to_tree(S_array::Array{SplitList{T},1}, t::Tree; strict=false, MCCs=nothing) where T
       out = SplitList(
               sort(collect(keys(t.lleaves))),
               Array{Split,1}(undef,0),
               ones(Bool, length(t.lleaves)), Dict{T, Split}(),
       )
       treesplits = SplitList(t)
       for S in S_array
               mS = map_splits_to_tree(S, t, treesplits; strict=strict, MCCs=MCCs)
               for s in mS
                       push!(out.splits, s)
               end
       end
       return out
end

"""
	map_splits_to_tree(S::SplitList, t::Tree)

Map splits `S` from another tree to `t`:
- restrain each split of `S` to `S.mask`
- find the corresponding internal node in `t`
- compute the split in `t` defined by this internal node.

Useful for resolving a tree with splits of another.
"""
function map_splits_to_tree(S::SplitList, t::Tree; strict=false, MCCs=nothing)
	treesplits = SplitList(t)
	return map_splits_to_tree(S, t, treesplits; strict=strict, MCCs=MCCs)
end
function map_splits_to_tree(S::SplitList, tree::Tree, treesplits::SplitList; strict=false, MCCs=nothing)
	mS = SplitList(
		S.leaves,
		Array{Split,1}(undef,0),
		ones(Bool, length(treesplits.leaves)),
		Dict{eltype(S.leaves), Split}()
	)
	for i in 1:length(S)
		ms = _map_split_to_tree(S, i, tree, treesplits; strict=strict, MCCs=MCCs)
		push!(mS.splits, ms)
	end
	return mS
end


#=
Map split `S[i]` to `t`.
=#
function _map_split_to_tree(S::SplitList, i::Integer, t::Tree, treesplits::SplitList; strict=false, MCCs=nothing)

	ms = Split(0)
    # Not lca in case lca(t, leaves(S,i)) contains extra leaves not in S[i]
    roots = TreeTools.blca([t.lleaves[x] for x in leaves(S,i)]...)
    if strict == true && !isnothing(MCCs)
        mcc_map = leaf_mcc_map(MCCs)
        assign_mccs!(mcc_map, t) 
        mcc_ = [r.data.dat["mcc"] for r in roots]
        mcc_ =  mcc_[mcc_.!=nothing]
        @assert !isnothing(mcc_)
        children = TreeTools.lca([t.lleaves[x] for x in leaves(S,i)]...).child
        sisters = children[children .∉ Ref(roots)]
        if any([(r.isroot || !(r.anc.data.dat["mcc"]==r.data.dat["mcc"] || r.data.dat["mcc"] == mcc_)) for r in sisters])
            return ms
        end
    end
    
	for r in roots
		if r.isleaf # `treesplits.splitmap` (probably) does not contain leaf-splits
			TreeTools.joinsplits!(ms, Split([findfirst(==(r.label), S.leaves)]))
		else
			TreeTools.joinsplits!(ms, treesplits.splitmap[r.label])
		end
	end
	return ms
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
