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
get_trees(no_trees, no_lineages; remove=false, debug=false, c=0.75, ρ = 0.05

returns trees as well as true MCCs
simulate a total of `no_trees` trees using the ARGTools package, 
specifying the `lineage_no` determines the number of nodes in each tree
remove - if internal branches should be removed (i.e. if trees should not be fully resolved, see parameter c)
c - Parameter to desribe how resolved trees are
ρ - Reassortment rate scaled to coalescence rate

"""
function get_trees(no_trees, no_lineages; remove=false, c=0.75, ρ = 0.05, N = 10_000)
    # Parameters of the ARG simulation
    N = N # pop size
    n = no_lineages # Number of lineages
    simtype = :flu
    r = get_r(ρ, n, N, simtype) # Absolute reassortment rate

    # Simulating the ARG
    arg = ARGTools.SimulateARG.simulate(N, r, n; K=no_trees, simtype);
    # The trees for the 2 segments
    trees = ARGTools.trees_from_ARG(arg; node_data = TreeTools.MiscData);
    if remove
        trees = remove_branches(trees; c=c, N = N)
    else
        trees = [convert(TreeTools.Tree{TreeTools.MiscData}, t) for t in trees]
    end
    return trees, arg
end

function remove_branches(input_trees; c=0.75, N = 10_000)
    trees = [copy(t) for t in input_trees]
    no_trees = length(trees)
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
        end
    end
    trees = [convert(TreeTools.Tree{TreeTools.MiscData}, t) for t in trees]
    return trees
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

function compute_consecutive_MCC_set(trees::Vector{TreeTools.Tree{T}}) where T
    M = TreeKnit.MCC_set(length(trees), [t.label for t in trees])
    for i in 1:(M.no_trees-1)
        for j in (i+1):M.no_trees
            oa = TreeKnit.OptArgs(;γ = 2., likelihood_sort = true, resolve = true,
                nMCMC = 25, verbose=false,)
            mCCs = TreeKnit.runopt(oa, trees[i], trees[j]; output = :mccs)
            TreeKnit.add!(M, mCCs, (i, j))
        end
    end
    return M
end
