
"""
get_mcc_map(MCCs::Vector{Vector{String}}, get_cluster_no =true)

Returns a dictionary of which MCC each leaf is in, MCCs are indexed by location in the Vector. 
"""
function leaf_mcc_map(MCCs::Vector{Vector{String}})
    mcc_map = Dict{String, Int}()
    for (i,mcc) in enumerate(MCCs)
        for node in mcc
            mcc_map[node] = i
        end
    end
    return mcc_map
end

"""
PRT!(n::TreeNode)

Assign `mcc`s to branches (i.e. their child node) by a Pre-order traversal starting at the root node `n`.
"""
function PRT!(n::TreeNode, k::Int)
	
    if isroot(n)
        n.data["mcc"] = []
        for pos in 1:k
            if !isempty(n.data["child_mccs"][pos]) && length(n.data["child_mccs"][pos])==1
                append!(n.data["mcc"], pop!(n.data["child_mccs"][pos]))
            else
                append!(n.data["mcc"], [nothing])
            end
        end
    else
        n.data["mcc"] =[]
        for pos in 1:k
            if n.anc.data["mcc"][pos] in n.data["child_mccs"][pos] # parent MCC part of children -> that is the MCC
                append!(n.data["mcc"], n.anc.data["mcc"][pos])
            elseif length(n.data["child_mccs"][pos])==1  # child is an MCC
                append!(n.data["mcc"], pop!(n.data["child_mccs"][pos]))
            else # no unique child MCC and no match with parent -> not part of an MCCs
                append!(n.data["mcc"], [nothing])
            end
        end
    end

	#delete!(n.data.dat, "child_mccs")

	if !isempty(n.child)
		for c in n.child
			PRT!(c, k)
		end
	end


end

function PRT!(t::Tree, k::Int)
	PRT!(t.root, k)
end

function PRT!(n::TreeNode)
    if isroot(n)
        if !isempty(n.data["child_mccs"]) && length(n.data["child_mccs"])==1
            n.data["mcc"] = pop!(n.data["child_mccs"])
        else
            n.data["mcc"] = nothing
        end
    else
        if n.anc.data["mcc"] in n.data["child_mccs"] # parent MCC part of children -> that is the MCC
            n.data["mcc"] = n.anc.data["mcc"]
        elseif length(n.data["child_mccs"])==1  # child is an MCC
            n.data["mcc"] = pop!(n.data["child_mccs"])
        else # no unique child MCC and no match with parent -> not part of an MCCs
            n.data["mcc"] = nothing
        end
    end
    #delete!(n.data.dat, "child_mccs")

    if !isempty(n.child)
        for c in n.child
            PRT!(c)
        end
    end

end

function PRT!(t::Tree)
	PRT!(t.root)
end

"""
	assign_mccs!(mcc_map::Dict{String, Int}, t::Vector{Tree{TreeTools.MiscData}})

Assign each node, (leaf and internal) node to the MCC that they are a part, takes dictionary `mcc_map` 
(leaf => MCC) as input and a tree `t`, if there is a conflict or it is unknown which MCC a node is part of 
this is labeled as None, this is a distinction to `mcc_map` which will throw an error if there is a conflict
between nodes, meaning that it cannot be called on unresolved trees. 

"""
function assign_mccs!(mcc_map::Dict{String, Int}, tree::Tree{TreeTools.MiscData}) 
	
	# assign MCCs to leaves
    for leaf in tree.lleaves
        leaf.second.data["child_mccs"] = Set([mcc_map[leaf.second.label]])
        leaf.second.data["mcc"] = mcc_map[leaf.second.label]
    end

    # reconstruct MCCs with Fitch algorithm
    for n in POT(tree)
        if !n.isleaf
            common_mccs = intersect([c.data["child_mccs"] for c in n.child]...)
            if !isempty(common_mccs)
                n.data["child_mccs"] = common_mccs
            else
                n.data["child_mccs"] = union([c.data["child_mccs"] for c in n.child]...)
            end
        end
    end

    PRT!(tree)
end
