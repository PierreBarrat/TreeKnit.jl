using ARGTools
using TreeTools

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
function get_trees(no_trees, no_lineages; get_real_MCCs=false)
    # Parameters of the ARG simulation
    N = 10_000 # pop size
    n = no_lineages # Number of lineages
    ρ = 0.1 # Reassortment rate scaled to coalescence rate
    simtype = :kingman
    r = get_r(ρ, n, N, simtype) # Absolute reassortment rate

    # Simulating the ARG
    arg = ARGTools.SimulateARG.simulate(N, r, n; K=no_trees, simtype);
    # The trees for the 2 segments
    trees = ARGTools.trees_from_ARG(arg; node_data = TreeTools.MiscData);

    if get_real_MCCs
        rMCCs = Vector{Vector{String}}()
        for k in 2:no_trees
            k_iters = Combinatorics.combinations(1:no_trees, k)
            for combination in k_iters
                append!(rMCCs, ARGTools.MCCs_from_arg(arg, combination...));
            end
        end
        return trees, rMCCs
    end
    return trees
end