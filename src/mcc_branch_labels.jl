
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
assign_mccs_PR!(n::TreeNode)

Assign `mcc`s to branches (i.e. their child node) by a Pre-order traversal starting at the root node `n`.
"""
function assign_mccs_PR!(n::TreeNode)
    if isroot(n)
        if length(n.data["child_mccs"])==1
            n.data["mcc"] = collect(n.data["child_mccs"])[1]
        else
            n.data["mcc"] = nothing
        end
    else
        if n.anc.data["mcc"] in n.data["child_mccs"] # parent MCC part of children -> that is the MCC
            n.data["mcc"] = n.anc.data["mcc"]
        elseif length(n.data["child_mccs"])==1  # child is an MCC
            n.data["mcc"] = collect(n.data["child_mccs"])[1]
        else # no unique child MCC and no match with parent -> not part of an MCCs
            n.data["mcc"] = nothing
        end
    end
    #delete!(n.data.dat, "child_mccs")

    if !isempty(n.child)
        for c in n.child
            assign_mccs_PR!(c)
        end
    end

end

function assign_mccs_PR!(t::Tree)
	assign_mccs_PR!(t.root)
end

"""
	assign_mccs!(t::Vector{Tree{TreeTools.MiscData}}, mcc_map::Dict{String, Int})

Assign each node, (leaf and internal) node to the MCC that they are a part, takes dictionary `mcc_map` 
(leaf => MCC) as input and a tree `t`, if there is a conflict or it is unknown which MCC a node is part of 
this is labeled as None, this is a distinction to `mcc_map` which will throw an error if there is a conflict
between nodes, meaning that it cannot be called on unresolved trees. 

"""
function assign_mccs!(tree::Tree{TreeTools.MiscData}, mcc_map::Dict{String, Int}) 
	
	# assign MCCs to leaves
    for leaf in tree.lleaves
        leaf.second.data["child_mccs"] = Set([mcc_map[leaf.second.label]])
        leaf.second.data["mcc"] = mcc_map[leaf.second.label]
    end

    # reconstruct MCCs with Fitch algorithm
    for n in POT(tree)
        if isroot(n)
            n.data["child_mccs"] = intersect([c.data["child_mccs"] for c in n.child]...)
        else
            if !n.isleaf
                common_mccs = intersect([c.data["child_mccs"] for c in n.child]...)
                if !isempty(common_mccs)
                    n.data["child_mccs"] = common_mccs
                else
                    n.data["child_mccs"] = union([c.data["child_mccs"] for c in n.child]...)
                end
            end
        end
    end

    assign_mccs_PR!(tree)
end
