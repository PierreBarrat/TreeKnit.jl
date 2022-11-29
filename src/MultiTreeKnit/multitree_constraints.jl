"""
	MCC_join_constraint(MCCs::Vector{Vector{Vector{String}}}; dict=false)

Join input lists of MCCs by recursively calculating their intersection.
"""
function MCC_join_constraint(MCCs::Vector{Vector{Vector{String}}}; dict=false)
    sets = [union([Set([m... ]) for m in mcc]...) for mcc in MCCs]
    @assert union(sets...) == sets[1] ## make sure labels are the same in all trees
    mcc_map = Dict{String, Vector{Int}}()
    reverse_mcc_map = Dict{Vector{Int}, Set{String}}()
    for MCC in MCCs
        for (i,mcc) in enumerate(MCC)
            for node in mcc
                if haskey(mcc_map, node)
                    append!(mcc_map[node], i)
                else
                    mcc_map[node] = [i]
                end
            end
        end
    end
    for node in keys(mcc_map)
        if haskey(reverse_mcc_map, mcc_map[node])
            push!(reverse_mcc_map[mcc_map[node]], node)
        else
            reverse_mcc_map[mcc_map[node]] = Set([node])
        end
    end
    if dict
        MCCs_new = Dict{Int, Set{String}}()
        for (num, key) in enumerate(keys(reverse_mcc_map))
            MCCs_new[num] = reverse_mcc_map[key]
        end
        return MCCs_new
    else
        MCCs_new = Vector{String}[]
        for (num, mcc) in reverse_mcc_map
            append!(MCCs_new, [sort(collect(mcc))])
        end
        return TreeKnit.sort(MCCs_new; lt=TreeKnit.clt)
    end
end

"""
	is_MCC_subset(MCC1::Vector{Vector{String}}, MCC2::Vector{Vector{String}})

function to check that every set in MCC1 is a subset of a set in MCC2
"""
function is_MCC_subset(MCC1::Vector{Vector{String}}, MCC2::Vector{Vector{String}})
    for mcc1 in MCC1
        subset = false
        for mcc2 in MCC2
            if issubset(Set{String}(mcc1), Set{String}(mcc2))
                subset = true
                break
            end
        end
        if !subset
            return false
        end
    end
    return true
end

function is_MCC_subset(MCC1::Dict{Int, Set{String}}, MCC2::Dict{String, Union{Int, Nothing}})
    for (keys, mcc1) in MCC1
        if length(Set([MCC2[m] for m in mcc1]))!=1
            return false
        end
    end
    return true
end


"""
    map_shared_branches!(MCC::Union{Nothing, Vector{Vector{String}}}, tree::Tree{TreeTools.MiscData})

Add a `shared_branch` parameter to all `trees`, branches with a `shared_branch` do not have a recombination event 
occuring on them as they connect clades that should be together according to the MCCs. If trees are not fully resolved 
and two mccs share a common ancestor this ancestor will be marked as `None` by the `assign_mccs!` function. If the 
`up_to_resolution` flag is used branches leading to such a shared node will still be marked as `shared_branch` 
(i.e. node.anc's `mcc` is marked as `None` but node has a sister with the same `mcc`).

The function proceeds by allocating each node to the MCC it should be in using the Fitch algorithm, 
then branches which are in a MCC with 2 or more nodes are marked with `shared_branch`.
"""
function map_shared_branches!(MCC::Union{Nothing, Vector{Vector{String}}}, tree::Tree{TreeTools.MiscData}; up_to_resolution=true)
	
	if isnothing(MCC)
		return nothing
	end

    shared_branches_map_ = map_shared_branches(MCC, tree; up_to_resolution)

	# assign MCCs to leaves
    for n in POT(tree)
        n.data["shared_branch"] = shared_branches_map_[n.label]
    end
    return nothing
end

function map_shared_branches!(MCC, tree; up_to_resolution=true)
@error "`map_shared_branches!` is only for trees with data attached to nodes, *i.e.* `Tree{TreeTools.MiscData}`.
	 Try `map_shared_branches` to get a map from nodes to shared branches.
	"
	error("Incorrect method")
end

"""
    map_shared_branches(MCC::Union{Nothing, Vector{Vector{String}}}, tree::Tree)

Create a dictionary which maps each node label to if that branch is shared (Bool) in the MCC. If trees are not fully
resolved and two mccs share a common ancestor this ancestor will be marked as `None` by the `assign_mccs` function. If the 
`up_to_resolution` flag is used branches leading to such a shared node are `shared`.

The function proceeds by allocating each node to the MCC it should be in using the Fitch algorithm, 
then branches which are in a MCC with 2 or more nodes are marked with `shared_branch`.
"""
function map_shared_branches(MCC::Union{Nothing, Vector{Vector{String}}}, tree::Tree; up_to_resolution=true)
    if isnothing(MCC)
		return nothing
	end
    shared_branches_map = Dict{String, Bool}()
	    mcc_map_ = map_mccs(tree, MCC)
        
    for n in POT(tree)
        if !isnothing(mcc_map_[n.label]) && length(MCC[mcc_map_[n.label]]) >1
            m = mcc_map_[n.label]
            if (isroot(n) || m==mcc_map_[n.anc.label])
                shared_branches_map[n.label] = true
            elseif up_to_resolution && (any([(c!=n && (mcc_map_[c.label]==m ||  any([mcc_map_[l.label] ==m for l in POTleaves(c)]))) for c ∈ n.anc.child]))
                shared_branches_map[n.label] = true
            else
                shared_branches_map[n.label] = false
            end
        else
            shared_branches_map[n.label] = false
        end
    end

    return shared_branches_map
end
