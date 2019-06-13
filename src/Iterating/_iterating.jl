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

and text files containing strains that have been removed, with two columns: strain name, and MCC name (so I can group them) 
=#

## A BIG problem that I see right now
# If I do not re-infer trees, the jointtree will quickly be completely irrelevant
# But I need it to adjust branch lengths! 
# I could re-infer only the joint tree, but i'm afraid that might change its MCCs... 

global segments = ["ha","na"]
global root_strain = "A/Victoria/361/2011"


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
"""
function prunetrees_mut(segtrees, jointtree ; verbose = true)

	## Resolving, finding MCCs, and adjusting branch length to that of the joint tree
	# resolving
	_resolve_trees!(segtrees, jointtree);
	# MCCs
	MCC = maximal_coherent_clades(collect(values(segtrees)))
	println("\n### MCCs ###\n")
	println("Found $(length(MCC)) MCCs of average size $(mean([length(x) for x in MCC]))")
	# Adjusting branch length in MCCs
	_adjust_branchlength!(segtrees, jointtree, MCC)

	## Computing ancestral states
	run_treetime(segtrees)

	## Counting mutations
	map(s->fasta2tree!(segtrees[s], "CrossMapping/ancestral_tree_$(s)_aln_$(s)/ancestral_sequences.fasta"), segments)
	segmd_ref = Dict(s=>make_mutdict!(segtrees[s]) for s in segments)
	crossmuts = compute_crossmuts(segtrees, MCC, segmd_ref)

	## Finding suspicious MCCs
	fmcc, tofilter, confidence = find_suspicious_mccs(crossmuts, segtrees)
	iroot = findall(x->x==root_strain, tofilter)
	for i in iroot
		splice!(tofilter, i)
		splice!(confidence, i)
	end
	return tofilter, fmcc, confidence
end



end