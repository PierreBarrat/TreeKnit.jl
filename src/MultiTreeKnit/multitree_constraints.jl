

"""
MCC_join_constraint(MCCs::Vector{Vector{Vector{String}}}; dict=false)

given as input a list of MCCs this function joins these MCC lists (sets) by calculating their intersection recursively
and returning their intersection
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
mark_shared_branches!(MCC::Union{Nothing, Vector{Vector{String}}}, trees::Vararg{Tree})

Add a `shared_branch` parameter to all `trees`, branches with a `shared_branch` do not have a recombination event 
occuring on them as they connect clades that should be together according to the MCCs. If trees are not fully resolved 
and two mccs share a common ancestor this ancestor will be marked as `None` by the `assign_mccs!` function. If the 
`up_to_resolution` flag is used branches leading to such a shared node will still be marked as `shared_branch` 
(i.e. node.anc's `mcc` is marked as `None` but node has a sister with the same `mcc`).

The function proceeds by allocating each node to the MCC it should be in using the Fitch algorithm, 
then branches which are in a MCC with 2 or more nodes are marked with `shared_branch`.
"""
function mark_shared_branches!(MCC::Union{Nothing, Vector{Vector{String}}}, trees::Vararg{Tree}; up_to_resolution=true)
	
	if isnothing(MCC)
		return nothing
	end

	# assign MCCs to leaves
    for tree in trees
	    map_mccs!(tree, MCC)

		for n in POT(tree)
			if !isnothing(n.data["mcc"]) && length(MCC[n.data["mcc"]]) >1
                m = n.data["mcc"]
                if (isroot(n) || m==n.anc.data["mcc"])
				    n.data["shared_branch"] = true
                elseif up_to_resolution && (any([(c!=n && (c.data["mcc"]==m ||  any([l.data["mcc"] ==m for l in POTleaves(c)]))) for c ∈ n.anc.child]))
                    n.data["shared_branch"] = true
                else
                    n.data["shared_branch"] = false
                end
            else
				n.data["shared_branch"] = false
			end
		end
	end
end
