module Iterating

using StatsBase, DelimitedFiles
using TreeTools, RecombTools
using RecombTools.SplitGraph


#= 
General structure of an iteration folder Itn 
3 directories
- InData
-- InData/msas
-- InData/trees
- AfterFilter
-- AfterFilter/msas
-- AfterFilter/trees
- CrossMapping

and text files containing strains that have been removed, with two columns: strain name, and MCC number (so I can group them) 
=#

## A BIG problem that I see right now
# If I do not re-infer trees, the jointtree will quickly be completely irrelevant
# But I need it to adjust branch lengths! 
# I could re-infer only the joint tree, but i'm afraid that might change its MCCs... 

# --> Reinfer a tree for each MCC, using the joint tree! 
# Use Interfacing.scalebranches for this

global segments = ["ha","na"]
global root_strain = "A/Victoria/361/2011"

let state = 0
	global mcc_idx() = (state += 1)
end

function run_it(it::Int64, cwd::String, segtrees_, jointtree_ ; mut=true, sa=true)
	segtrees = deepcopy(segtrees_)
	jointtree = deepcopy(jointtree_)

	## Arrange directories
	make_directories(it, cwd)

	## Mutation based pruning
	if mut
		torm_strains_mut, torm_mccs_mut = ... 
	else
		torm_strains_mut = Array{String,1}(undef,0)
		torm_mccs_mut = Array{String,1}(undef,0)
	end

	for (s,t) in segtrees
		segtrees[s] = prunenodes(t, torm_strains_mut)
	end
	jointtree = prunenodes(jointtree, torm_strains_mut)

	## Topology based pruning
	if sa
		torm_strains_top, torm_mccs_top = ...
	else
		torm_strains_sa = Array{String,1}(undef,0)
		torm_mccs_sa = Array{String,1}(undef,0)
	end

	for (s,t) in segtrees
		segtrees[s] = prunenodes(t, torm_strains_sa)
	end
	jointtree = prunenodes(jointtree, torm_strains_sa)

	## Export (also filters fasta alignments)
	write_results(it, cwd, segtrees, jointtree, torm_strains_mut, torm_strains_sa) 

	return segtrees, jointtree
end

"""
	removeMCCs!(segtrees, mcclist)

Remove all strains in `mcclist` from `segtrees`. `mcclist` should be an array of arrays of strings. (*i.e.* one element for each MCC to be removed). 
Removed strains are written to the end of `outfile`. 
"""
function removeMCCs!(segtrees, mcclist ; outfile = "StrainsToFilter.txt")
	out = 
	for m in mcclist
	end
end


"""
	_resolve_trees!(segtrees)

Delete null branches above internal nodes. Then, pairwise resolving of trees in `segtrees`. 
Since I am not re-inferring trees, I cannot use the joint tree to resolve. 


### Notes
Initial null (*i.e.* too short to bear a mutation) branches above internal nodes are not trusted. --> They are set to 0, corresponding internal nodes are removed.  
Trees are resolved using the joint tree, introducing new very small branches. These are trusted because of topological evidence.  
Finally, leaves that have too short of a branch are also stretched to a minimum threshold.  
"""
function _resolve_trees!(segtrees ; verbose=false)
	# Removing null branches
	for (s,t) in segtrees
		L = length(t.leaves[1].data.sequence)
		delete_null_branches!(t.root, threshold = 1/L/3)
		t = node2tree(t.root)
		t = remove_internal_singletons(t)
		segtrees[s] = t
	end

	# Pairwise resolving
	for (s,t) in segtrees
		for (s2,t2) in segtrees
			if s2 != s
				t = resolve_trees(t, t2, rtau = 1/L/3, verbose=verbose)
			end
		end
		resolve_null_branches!(t, tau = 1/L/3)
		segtrees[s] = t
	end
end





include("mutbased.jl")
include("topologybased.jl")

end