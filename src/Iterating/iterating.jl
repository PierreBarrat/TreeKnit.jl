module Iterating

using StatsBase
using DelimitedFiles
using TreeTools, RecombTools

global segments = ["ha","na","pb1","pb2","pa"]
global root_strain = "A/Victoria/361/2011"

"""
"""
function run_it(it::Int64)
	# Setting up directories
	println("\n\n ### ITERATION $it ###\n")
	mkpath("It$(it)")
	cd("It$(it)")
	mkpath("InData")
	mkpath("OutData")
	# Copy data from previous It. 
	if it != 0
		for s in segments
			run(`cp ../It$(it-1)/OutData/outmsas/filtered_h3n2_$(s).fasta InData/h3n2_$(s).fasta`)
		end
	end
	cd("..")
	# Running Augur
	runaugur(it)
	# Filtering strains
	prunetrees_sup(it)
end

"""
"""
function runaugur(it::Int64)
	# Setting directories
	run(`cp -r AugurData_template/ It$(it)/AugurData`)
	map(x->run(`cp -t It$(it)/AugurData/results It$(it)/InData/h3n2_$(x).fasta`), segments)

	# Running the augur script
	cd("It$(it)/AugurData")
	run(`bash pierres_scripts/data_processing`)
	run(`julia pierres_scripts/concatenate.jl `)
	run(`bash pierres_scripts/process_allseg`)
	# Copying to OutData
	cd("../OutData")
	mkpath("msas")
	mkpath("trees")
	for s in (segments..., "all")
		run(`cp ../AugurData/results/aligned_simple_h3n2_$(s).fasta msas`)
		run(`cp ../AugurData/results/tree_refined_h3n2_$(s).nwk trees`)
	end
	cd("../../")
end

"""
- Resolve trees
- Find MCCs
- Compute ancestral states
- Find suspicious MCCs
- Filter sequences
"""
function prunetrees_mut(it::Int64)
	cd("It$(it)/OutData")
	## Reading trees
	segtrees = read_trees(segments)
	jointtree = read_trees(("all", ))["all"];

	## Resolving
	_resolve_trees!(segtrees, jointtree);


	## MCCs
	MCC = maximal_coherent_clades(collect(values(segtrees)))
	println("\n### MCCs ###\n")
	println("Found $(length(MCC)) MCCs of average size $(mean([length(x) for x in MCC]))")

	## Adjusting branch length in MCCs
	_adjust_branchlength!(segtrees, jointtree, MCC)

	## Writing resolved trees for auspice
	# This is not very well placed in the overall module ... move it? 
	write_trees(segtrees)
	write_trees(Dict("all"=>jointtree))

	## I SHOULD RE-EXPORT FOR AUSPICE HERE

	## Computing ancestral states
	println("\n### Running treetime...###")
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
	println("\n### FILTERING ###")
	println("Filtering $(length(tofilter)) strains ($(length(fmcc)) MCCs).")
	writedlm("StrainsToFilter_mut.txt",cat(tofilter, confidence, dims=2))
	mkpath("outmsas")
	for s in segments
		filter_fasta("../AugurData/results/aligned_h3n2_$(s).fasta", "outmsas/filtered_h3n2_$(s).fasta", tofilter)
	end

	# 
	cd("../../")
end

"""
- Find MCCs
- Find supra-MCC `s` for each MCC (initializing with `first(segtrees)`)
- Attempt to remove single MCCs from `s`
- Attempt to remove pairs of MCCs from `s`
- If one of the two above procedures completely resolved `s`, store the corresponding removed MCCs in a `toremove` list. 
- Filter MCCs in `toremove` from the dataset. 
"""
function prunetrees_sup(it::Int64)
	cd("It$(it)/OutData")
	## Reading trees
	segtrees = read_trees(segments)
	jointtree = read_trees(("all", ))["all"];

	## Resolving
	_resolve_trees!(segtrees, jointtree);


	## MCCs
	MCC = maximal_coherent_clades(collect(values(segtrees)))
	println("\n### MCCs ###\n")
	println("Found $(length(MCC)) MCCs of average size $(mean([length(x) for x in MCC]))")

	## Adjusting branch length in MCCs
	_adjust_branchlength!(segtrees, jointtree, MCC)

	## Writing resolved trees for auspice 
	# This is not very well placed in the overall module ... move it? 
	write_trees(segtrees)
	write_trees(Dict("all"=>jointtree))

	## I SHOULD RE-EXPORT FOR AUSPICE HERE

	##
	supra = RecombTools.supraMCCs(values(segtrees), MCC);
	mcc_scores_1, mcc_scores_2 = _compute_mcc_scores(segtrees, jointtree, supra)
	##

	## Finding suspicious MCCs
	fmcc, tofilter = find_suspicious_mccs_sup(segtrees, mcc_scores_1, mcc_scores_2)
	iroot = findall(x->x==root_strain, tofilter)
	for i in iroot
		splice!(tofilter, i)
	end

	# Filtering
	println("\n### FILTERING ###")
	println("Filtering $(length(tofilter)) strains ($(length(fmcc)) MCCs).")
	writedlm("StrainsToFilter_supra.txt",tofilter)
	mkpath("outmsas")
	for s in segments
		filter_fasta("../AugurData/results/aligned_h3n2_$(s).fasta", "outmsas/filtered_h3n2_$(s).fasta", tofilter)
	end

	# 
	cd("../../")
end

"""
Read trees from file
"""
function read_trees(segments)
	out = Dict(s=>read_tree("trees/tree_refined_h3n2_$s.nwk") for s in segments)
	map(s->fasta2tree!(out[s], "msas/aligned_simple_h3n2_$(s).fasta"), segments)
	return out
end

"""
Delete null branches above internal nodes. Then resolve all trees in `segtrees` using `jointtree`. 
### Notes
Initial null (*i.e.* too short to bear a mutation) branches above internal nodes are not trusted. --> They are set to 0, corresponding internal nodes are removed.  
Trees are resolved using the joint tree, introducing new very small branches. These are trusted because of topological evidence.  
Finally, leaves that have too short of a branch are also stretched to a minimum threshold.  
"""
function _resolve_trees!(segtrees, jointtree; verbose=false)
	verbose && println("\n### RESOLVING ###\n")
	delete_null_branches!(jointtree.root, threshold = 1/length(jointtree.leaves[1].data.sequence)/10)
	jointtree = node2tree(remove_internal_singletons(jointtree).root)
	check_tree(jointtree)
	# Joint tree as a reference
	for (k,t) in segtrees
		verbose && println("$k...")
		L = length(t.leaves[1].data.sequence)
		delete_null_branches!(t.root, threshold = 1/L/10)
		t = node2tree(t.root)
		t = resolve_trees(t, jointtree, rtau = 1/L/10, verbose=verbose)
		t = remove_internal_singletons(t)
		resolve_null_branches!(t, tau = 1/L/10)
		segtrees[k] = t
	end
end

"""
Write tree, as well as  branch length files for viewing with auspice. 
"""
function write_trees(tdict)
	for (k,t) in tdict
		write_newick("trees/tree_resolved_h3n2_$(k).nwk", t)
		write_branchlength("trees/branchlengths_resolved_$(k).json", t, "results/aligned_simple_h3n2_$(k).fasta", "../trees/tree_resolved_h3n2_$(k).nwk")
	end

end

"""
"""
function _adjust_branchlength!(segtrees, jointtree, MCC)
	name_mcc_clades!((collect(values(segtrees))..., jointtree), MCC)
	adjust_branchlength!(collect(values(segtrees)), jointtree, MCC)
end

"""
"""
function run_treetime(segtrees)
	mkpath("CrossMapping")
	cd("CrossMapping")
	run(`cp ../../../treetime_script_template.sh treetime_script`)
	for (s,t) in segtrees
		write_newick("tree_$(s)_adjusted.nwk", t)
	end
	map(k->run(`cp ../msas/aligned_simple_h3n2_$(k).fasta ./`), collect(keys(segtrees)))
	run(`bash treetime_script`)
	cd("..")
end

"""
Construct a dictionary of the form `"mcc label"=>CrossMutations`. 
"""
function compute_crossmuts(segtrees, MCC, segmd_ref)
	# Building a [tree, aln] dictionary of mutations
	crossmuts_all = Dict()
	for a in keys(segtrees)
	    for t in keys(segtrees)
	        crossmuts_all[t,a] = parse_nexus("CrossMapping/ancestral_tree_$(t)_aln_$(a)/annotated_tree.nexus")
	    end
	end

	# Inverting it to have a "per MCC" dictionary of mutations
	crossmuts_mccs = Dict()
	for m in MCC
	    lab = lca([segtrees[segments[1]].lnodes[x] for x in m]).label
	    cm = CrossMutations()
	    for a in keys(segtrees)
	        for t in keys(segtrees)
	            cm.crossmut[t,a] = crossmuts_all[t,a][lab]
	            cm.suspicious[t,a] = 0
	            for mut in cm.crossmut[t,a]
	                if get(segmd_ref[a][1],mut,0)==0 && mut[2]!=5 && mut[3]!=5
	                    cm.suspicious[t,a] += 1
	                end
	            end
	        end
	    end
	    crossmuts_mccs[lab] = cm
	end	
	return crossmuts_mccs
end

"""
"""
function find_suspicious_mccs(crossmuts, segtrees)
	fmcc = []
	tofilter = []
	confidence = []
	for (l,cm) in crossmuts
	    if sum(values(cm.suspicious))>0
	        push!(fmcc, l)
	        push!(tofilter, [x.label for x in node_leavesclade(segtrees[segments[1]].lnodes[l])]...)
        	push!(confidence, repeat([sum(collect(values(cm.suspicious)))], length(node_leavesclade(segtrees[segments[1]].lnodes[l])))...)
	    end
	end
	return fmcc, tofilter, confidence
end

"""
"""
function filter_fasta(ifasta, ofasta, filt)
	f = readlines(ifasta)
	open(ofasta, "w") do of
		flag = true
		for (i,l) in enumerate(f)
			if l[1] == '>'
				flag = true
				for m in filt
					if match(Regex(m), l) != nothing
						flag = false
						break
					end
				end
				if flag
					write(of, "$l\n")
				end
			else
				if flag
					write(of, "$l\n")
				end
			end
		end
	end
end

"""
For individual MCCs, should return an array `score_ind` of length `length(supra)`. `score_ind[i]` is a dictionary with MCC labels as keys and score as value. 
"""
function _compute_mcc_scores(segtrees, jointtree, supra)
	score_ind = []
	score_pairs = []
	for (i,sup) in enumerate(supra)
	    suptrees = Dict{String, Tree}()
	    sup_jointtree = node2tree(prunenode(lca([jointtree.lnodes[x] for x in sup.labels])))
	    for (s,t) in segtrees
	        suptrees[s] = node2tree(prunenode(lca([t.lnodes[x] for x in sup.labels])))
	    end
	    MCCloc = maximal_coherent_clades(collect(values(suptrees)))
	    s1 = compute_mcc_scores(suptrees, sup_jointtree, MCCloc)
	    s2 = compute_mcc_scores_pairs(suptrees, sup_jointtree, MCCloc)
	    push!(score_ind, s1)
	    push!(score_pairs, s2)
	end
	return score_ind, score_pairs
end

"""
- s1: individual scores for given supra -- s2: pairwise scores for ... Those are dictionaries `s[m] = score`. The relevant score is `score[2]`
- if `s1[m][2] == 1` for some MCC `m`, `push!(toremove, m)`. break
"""
function find_suspicious_mccs_sup(segtrees, mcc_scores_1, mcc_scores_2)
	toremove = Set()
	for (i,s1) in enumerate(mcc_scores_1)
		s2 = mcc_scores_2[i]
		if !isempty(s1) && findmin([x[2] for x in values(s1)])[1] == 1 # there a unique an MCC which blocks resolving 
			map(x->push!(toremove, x), findall(x->x==1, Dict(k=>v[2] for (k,v) in s1)) )
		elseif !isempty(s2) && findmin([x[2] for x in values(s2)])[1] == 1
        	map(x->(push!(toremove,x[1]);push!(toremove,x[2])), findall(x->x==1, Dict(k=>v[2] for (k,v) in s2)))
    	end
    end
    tofilter = []
    for m in toremove
    	append!(tofilter, [x.label for x in node_leavesclade(segtrees["ha"].lnodes[m])])
    end
    return toremove, tofilter
end



end