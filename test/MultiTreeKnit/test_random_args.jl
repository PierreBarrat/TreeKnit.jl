using ARGTools

println("######## MultiTreeKnit on random ARGs ##########")

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
c - Parameter to describe how resolved trees are
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

@testset "check is_degenerate(MCCs) == (consistency_rate!=0.0) holds for random multitrees" begin
    repeat = 0
    while repeat < 10
        trees, arg = get_trees(3, 50, remove=true, c=0.75);
        label!(trees[1], "a")
        label!(trees[2], "b")
        label!(trees[3], "c")
        iMCCs = TreeKnit.get_infered_MCC_pairs!(trees, consistent = true, constraint_cost=4.)
        @test TreeKnit.is_degenerate(iMCCs) == (TreeKnit.consistency_rate(iMCCs, trees)!=0.0)
        repeat += 1
    end
end

function check_parallelized_TK()
    for i in 1:100
        true_trees, arg = get_trees(8, 15; ρ=(10^-0.1))
        labels = [t.label for t in true_trees]
        rMCCs = TreeKnit.MCC_set(8, labels, TreeKnit.get_real_MCCs(8, arg))
        for no_trees in 2:8
            rand_order = sample(1:8, no_trees, replace = false)
            unresolved_trees = [copy(t) for t in true_trees[rand_order]]
            unresolved_trees = remove_branches(unresolved_trees; c=0.75)

            i_trees = [copy(t) for t in unresolved_trees]
            try
                fc_i_MCCs = TreeKnit.get_infered_MCC_pairs!(i_trees, TreeKnit.OptArgs(consistent = true, parallel=true))
            catch e
                print("n: "*string(no_trees))
                for t in true_trees
                    print_tree_ascii(" ", t)
                end
                throw(error())
            end
        end
    end
    return true
end

@testset "check parallelized TK runs correctly for random multitrees" begin
    @test check_parallelized_TK()
end