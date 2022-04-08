using ARGTools

using PyCall
using Conda
Conda.add("biopython")

## use a python wrapper to plot trees using Bio.Phylo
py"""
from Bio import Phylo
from io import StringIO
import matplotlib.pyplot as plt

def print_tree(tree_string, number):
    ax = plt.axes()
    tree = Phylo.read(StringIO(tree_string), "newick")
    Phylo.draw_ascii(tree)
    #plt.title("Tree "+ str(number))
    #Phylo.draw(tree, axes=ax)
"""

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
    get_trees(no_trees, no_lineages)
    
    simulate a total of `no_trees` trees using the ARGTools package, 
    specifying the `lineage_no` determines the number of nodes in each tree
"""
function get_trees(no_trees, no_lineages)
    # Parameters of the ARG simulation
    N = 10_000 # pop size
    n = no_lineages # Number of lineages
    ρ = 0.1 # Reassortment rate scaled to coalescence rate
    simtype = :kingman
    r = get_r(ρ, n, N, simtype) # Absolute reassortment rate

    tree_pairs = ceil(Int64, no_trees/2)
    last_tree_needed = (no_trees% 2 == 0)
    print(last_tree_needed)
    tree_strings = String[]
    trees = Any[]

    for i in 1:tree_pairs
        println(i)
        # Simulating the ARG
        arg = ARGTools.SimulateARG.simulate(N, r, n; simtype);
        # The trees for the 2 segments
        t_a, t_b = ARGTools.trees_from_ARG(arg; node_data = TreeTools.MiscData);

        tree_string= "";
        tree_string = TreeTools.write_newick!(tree_string, t_a.root)
        push!(tree_strings, tree_string)
        push!(trees, t_a)
        py"print_tree"(tree_string, 2*(i-1) +1)
        if i==tree_pairs && !last_tree_needed
            println("break")
            break
        end
        tree_string= "";
        tree_string = TreeTools.write_newick!(tree_string, t_b.root)
        push!(tree_strings, tree_string)
        push!(trees, t_b)
        py"print_tree"(tree_string, 2*(i-1) +2)
    end

    return trees, tree_strings
end