

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
