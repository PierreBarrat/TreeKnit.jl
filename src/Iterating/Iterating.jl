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
# global root_strain = "A/Victoria/361/2011"
global root_strain = "root"

let state = 0
	global mcc_idx() = (state += 1)
	global reset_mcc_idx() = (state = 0)
end
#=
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
=#

function run_it(it::Int64, segtrees_, jointMSA ; mut=true, sa=true, writetrees=true)
	println("############# It$it #############")
	segtrees = deepcopy(segtrees_)
	if it == 1
		reset_mcc_idx()
	end
	## Arrange directories
	mkpath("OutData/It$it")
	## Mutation based pruning
	if mut
		println("Pruning based on mutations")
		# Checking that sequences are present
		if !prod(prod(!isempty(x.data.sequence) for x in collect(values(segtrees[s].leaves))) for s in keys(segtrees))
			error("Trees do not have sequences attached")
		end
		segtrees, torm_mccs, mcc_names = prunetrees_mut(segtrees, jointMSA, cwd="OutData/It$it", verbose=true)
		mcclist = [mcc_names[x] for x in torm_mccs]
		filtered_strains_mut = removeMCCs!(segtrees, mcclist, outfile="OutData/It$it/FilteredStrains_mut.txt")
		println("Filtering $(length(mcclist)) MCCs, corresponding to $(length(filtered_strains_mut[:,1])) strains.")
	elseif sa
		println("Pruning based on topology")
		segtrees, mcclist = prunetrees_sa(segtrees, verbose=true, maxsize=75)
		filtered_strains_sa = removeMCCs!(segtrees, mcclist, outfile="OutData/It$it/FilteredStrains_sa.txt")
		println("Filtering $(length(mcclist)) MCCs, corresponding to $(length(filtered_strains_sa[:,1])) strains.")
	end
	#
	for (s,t) in segtrees
		write_newick("OutData/It$it/tree_$(s)_pruned.nwk", t.root)
	end
	#
	if mut && sa
		filtered_strains = [filtered_strains_mut, filtered_strains_sa]
	elseif sa
		filtered_strains = filtered_strains_sa
	elseif mut
		filtered_strains = filtered_strains_mut
	else
		filtered_strains = []
	end
	#
	return segtrees, filtered_strains
end


"""
	removeMCCs!(segtrees, mcclist)

Remove all strains in `mcclist` from `segtrees`. `mcclist` should be an array of arrays of strings. (*i.e.* one element for each MCC to be removed). 
Removed strains are written to the end of `outfile`. 
"""
function removeMCCs!(segtrees, mcclist ; outfile = "StrainsToFilter.txt")
	filtered = Array{Any, 2}(undef, 0, 2)
	if mapreduce(x->length(x), +, mcclist, init=0) == length(first(segtrees)[2].leaves)
		filtered = cat(collect(keys(first(segtrees)[2].lleaves)), repeat([mcc_idx()], length(first(segtrees)[2].leaves)),dims=2)
		for (s,t) in segtrees
			segtrees[s] = node2tree(TreeNode())
		end
	else
		for (i,m) in enumerate(mcclist)
			idx = mcc_idx()
			for n in m
				filtered = [filtered ; [n idx]]
			end
			for (s,t) in segtrees
				segtrees[s] = prunenode(t, m)
			end
			# println(i)
			# println(share_labels(segtrees["ha"], segtrees["na"]))
		end
		for (s,t) in segtrees
			segtrees[s] = remove_internal_singletons(t)
		end
	end
	writedlm(outfile, filtered)
	return filtered
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
		if L == 0 
			L = 1e4
		end
		delete_null_branches!(t.root, threshold = 1/L/3)
		t = node2tree(t.root)
		t = remove_internal_singletons(t)
		segtrees[s] = t
	end

	# Pairwise resolving
	for (s,t) in segtrees
		L = length(t.leaves[1].data.sequence)
		if L == 0
			L = 1e4
		end
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