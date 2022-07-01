using ARGTools
using TreeTools
using Random

# get recombination rate
function get_r(ρ, n, N, simtype::Symbol)
    if simtype == :kingman
        return ρ * n / N
    elseif simtype == :yule
        return ρ / N
    elseif simtype == :flu
    	return ρ * n^0.2 / N
    else
        @error "Unrecognized `simtype`."
    end
end

"""
    get_trees(no_trees, no_lineages, get_real_MCCs)
    
    returns trees as well as true MCCs
    simulate a total of `no_trees` trees using the ARGTools package, 
    specifying the `lineage_no` determines the number of nodes in each tree
"""
function get_trees(no_trees, no_lineages; remove=false, debug=false)
    # Parameters of the ARG simulation
    N = 10_000 # pop size
    n = no_lineages # Number of lineages
    ρ = 0.05 # Reassortment rate scaled to coalescence rate
    simtype = :kingman
    r = get_r(ρ, n, N, simtype) # Absolute reassortment rate

    # Simulating the ARG
    arg = ARGTools.SimulateARG.simulate(N, r, n; K=no_trees, simtype);
    # The trees for the 2 segments
    trees = ARGTools.trees_from_ARG(arg; node_data = TreeTools.MiscData);
    if remove
        c = 0.75
        for i in range(1, no_trees)
            tree = trees[i]
            delete_list = String[]
            for node in internals(tree)
                if !node.isroot
                    Pr = exp(-node.tau/(c*N))
                    if rand() <= Pr
                        push!(delete_list, node.label)
                    end
                end
            end
            for node_label in delete_list
                delete_node!(trees[i], node_label, ptau=true)
                if debug
                    println("removing node: "* node_label)
                end
            end
        end
    end
    return trees, arg
end

function get_real_MCCs(no_trees, arg)
    rMCCs = Vector{Vector{String}}[]
    for k in 2:no_trees
        k_iters = Combinatorics.combinations(1:no_trees, k)
        for combination in k_iters
            push!(rMCCs, ARGTools.MCCs_from_arg(arg, combination...));
        end
    end
    return rMCCs
end

function get_tree_order(trees ;order="RF_distance")
    no_trees = length(trees)
    if order=="resolution" ##start with most resolved trees
        resol_index = [resolution_value(t) for t in trees]
        permvec = sortperm(resol_index, rev=true)
    else
        ##start with trees that are most similar to all other trees
            RF_index = Float16[]
            for i in range(1, no_trees)
                rf = 0
                for j in range(1, no_trees)
                    if j!=i
                        rf += RF_distance(trees[i], trees[j])^2
                    end
                end
                push!(RF_index, rf/(no_trees-1))
            end
            permvec = sortperm(RF_index)
    end

    return permvec
end