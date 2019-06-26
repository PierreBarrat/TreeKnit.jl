module Iterating

using StatsBase, DelimitedFiles
using TreeTools, RecombTools
using RecombTools.SplitGraph

# global segments = ["ha","na","pb1","pb2","pa"]
global segments = ["ha","na"]
global root_strain = "A/Victoria/361/2011"


export read_trees, _resolve_trees!, _adjust_branchlength!, write_trees, prunetrees_mut, prunetrees_irem


"""
"""
function run_it(it::Int64 ; mut=true, sa=true, irem=false)
	# Setting up directories
	println("\n\n ### ITERATION $it ###\n")
	mkpath("It$(it)")
	cd("It$(it)")
	mkpath("InData")
	mkpath("OutData")
	# Copy data from previous It. 
	if it != 0
		println(pwd())
		for s in segments
			run(`cp ../It$(it-1)/OutData/outmsas/filtered_h3n2_$(s).fasta InData/h3n2_$(s).fasta`)
		end
	end
	cd("..")
	# Running Augur
	runaugur(it)
	# Filtering strains
	return prunetrees(it, mut=mut, sa=sa, irem=irem)
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
"""
function run_auspice_export(it::Int64)
	run(`bash pierres_scripts/set_branchlength_auspice.sh`)
end

"""
- Resolve trees
- Find MCCs
- Compute ancestral states
- Find suspicious MCCs
- Filter sequences
"""
function prunetrees(it::Int64 ; mut=true, sa=true,irem=false)
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
	cd("../AugurData")
	run_auspice_export(it)
	cd("../OutData")

	# Filtering
	tofilter_strains_mut = []
	tofilter_mccs_mut =[]
	tofilter_strains_sa = []
	tofilter_mccs_sa = []
	tofilter_strains_sup = []
	tofilter_mccs_sup = []
	if mut
		println("\n### Mutation based procedure...###")
		tofilter_strains_mut, tofilter_mccs_mut, confidence = prunetrees_mut(segtrees, MCC)		
	end
	if irem
		println("\n### Computing supra MCCs... ###")
		tofilter_strains_sup, tofilter_mccs_sup = prunetrees_irem(segtrees, jointtree, MCC)
	end
	if sa
		println("\n### Local optimization on supra MCCs... ###")
		# tofilter_strains_sa, tofilter_mccs_sa = prunetrees_sa(segtrees, maximal_coherent_clades(collect(values(segtrees))))
		tofilter_strains_sa, tofilter_mccs_sa = prunetrees_sa(segtrees, jointtree)
	end

	tofilter_strains_all = cat(tofilter_strains_mut,tofilter_strains_sup,tofilter_strains_sa,dims=1)
	tofilter_mccs_all = cat(tofilter_mccs_mut,collect(tofilter_mccs_sup),tofilter_mccs_sa,dims=1)

	println("\n### FILTERING ###")
	mut && println("Based on mutations: Filtering $(length(tofilter_strains_mut)) strains ($(length(tofilter_mccs_mut)) MCCs).")
	irem && println("Based on supra MCCs: Filtering $(length(tofilter_strains_sup)) strains ($(length(tofilter_mccs_sup)) MCCs).")
	sa && println("Based on local optimization: Filtering $(length(tofilter_strains_sa)) strains ($(length(tofilter_mccs_sa)) MCCs).")
	println("--> In total: Filtering $(length(tofilter_strains_all)) strains ($(length(tofilter_mccs_all)) MCCs).")

	filter_strains(tofilter_strains_all, tofilter_mccs_all)
	cd("../../")
	return tofilter_strains_all

end




"""
"""
function prunetrees_mut(segtrees, MCC)
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

"""
"""
function prunetrees_irem(segtrees, jointtree, MCC)
	supra, scores_1, scores_2, scores_3 = _compute_mcc_scores_irem(segtrees, jointtree, MCC)
	smcc = find_suspicious_mccs_irem(supra, scores_1, scores_2, scores_3)
	toremove = Array{String,1}(undef, 0)
	for m in smcc
		append!(toremove, node_leavesclade_labels(segtrees[first(segments)].lnodes[m]))
	end
	return toremove, smcc
end

"""
"""
function prunetrees_sa(segtrees, jointtree ; maxsize = 50, verbose=true, NIT = 100)
	st = deepcopy(segtrees)
	jt = deepcopy(jointtree)

	torm_strains = []
	torm_mccs = []
	# Init
	_resolve_trees!(st, jt);
	MCC = maximal_coherent_clades(collect(values(st)))
	_adjust_branchlength!(st, jt, MCC);

	Eold = SplitGraph.count_mismatches(collect(values(st))...)
	MCC_old = MCC
	nleaves_old = length(st["ha"].leaves)
	# Iterating
	for it in 1:NIT
	    println("-- SA it.$it --")
	    supra = unique([Set(x) for x in RecombTools.supraMCC(values(st), MCC)])
	    # Finding a good supra to work on
	    flag = false
	    for s in supra
	        subtrees = Dict{String, Tree}()
	        for (k,t) in st
	            r = lca([t.lnodes[x] for x in s])
	            if !r.isroot
	                subtrees[k] = node2tree(prunenode(r)[1])
	            else
	                subtrees[k] = deepcopy(st[k])
	            end
	        end
	        locMCC = maximal_coherent_clades(collect(values(subtrees)))
	        
	        # 
	        if length(locMCC) < 50
	            out = SplitGraph.opttrees(collect(values(subtrees))..., γ=0.95, Trange = .5:-0.01:0.05, M=2500);
	            torm = unique(cat(cat(out[1]..., dims=1)...,dims=1))
	            for s in Iterating.segments
	                st[s] = prunenodes(st[s], torm)
	            end
	            jt = prunenodes(jt, torm);
	
	            _resolve_trees!(st, jt,verbose=false);
	            MCC = maximal_coherent_clades(collect(values(st)))
	            if verbose
	            	println("Removing $(length(torm)) strains")
	            	println("Number of mismatches (old/new): ", Eold, " ", SplitGraph.count_mismatches(collect(values(st))...))
	            	println("Mean size (old/new): ", mean(length(x) for x in MCC_old), " ", mean(length(x) for x in MCC))
	            	println("# of MCCs (old/new): ", length(MCC_old), " ", length(MCC))
	            	println("# of strains (old/new): $nleaves_old $(length(st["ha"].leaves))")
	            end
	            append!(torm_strains, torm)
	            append!(torm_mccs, unique(cat(out[2]..., dims=1)))
	            flag = true
	            break
	        end
	    end
	    if !flag
	    	break
	    end
	end

	return torm_strains, torm_mccs
end

# """
# """
# function prunetrees_sa(segtrees, MCC ; maxsize=50, verbose=true)
# 	supra = unique(RecombTools.supraMCC(values(segtrees), MCC))
# 	torm_strains = Array{String,1}(undef,0)
# 	torm_mcc = Array{String,1}(undef,0)
# 	i=1
# 	c = 0
# 	for s in supra
# 		verbose && print("$i/$(length(supra))      	\r")
# 		i+=1
# 		# Making subtrees
# 		subtrees = Dict{String, Tree}()
# 		for (k,t) in segtrees
# 			r = lca([t.lnodes[x] for x in s])
# 			if !r.isroot
# 				subtrees[k] = node2tree(prunenode(r)[1])
# 			else
# 				subtrees[k] = deepcopy(segtrees[k])
# 			end
# 		end
# 		locMCC = maximal_coherent_clades(collect(values(subtrees))) # stupid way of checking for length. I should do something smarter

# 		if length(locMCC) < maxsize
# 			c += 1
# 			# Optimization
# 			out = SplitGraph.opttrees(collect(values(subtrees))..., γ=0.99, Trange = .5:-0.01:0.05, M=2500)
			
# 			if out[end]
# 				append!(torm_strains, cat(out[1]...,dims=1))
# 				append!(torm_mcc, out[2])
# 			else
# 				@warn "SA did not converge"
# 			end
# 		end
# 	end
# 	println("Found $c common clades amenable for local optimization")
# 	return unique(torm_strains), unique(torm_mcc)
# end

"""
- Find MCCs
- Find supra-MCC `s` for each MCC (initializing with `first(segtrees)`)
- Attempt to remove single MCCs from `s`
- Attempt to remove pairs of MCCs from `s`
- If one of the two above procedures completely resolved `s`, store the corresponding removed MCCs in a `toremove` list. 
- Filter MCCs in `toremove` from the dataset. 
"""
function _compute_mcc_scores_irem(segtrees, jointtree, MCC)
	scores_1 = []
	scores_2 = []
	scores_3 = []
	threshold_indiv = 50
	threshold_pair = 25
	threshold_triplets = 9

	# Reducing input trees to MCCs
	segtrees_reduced = Dict{String, Tree}()
	for (s,t) in segtrees
		segtrees_reduced[s] = reduce_to_mcc(t, MCC)
	end
	jointtree_reduced = reduce_to_mcc(jointtree, MCC)

	#
	supra = supraMCC(cat(collect(values(segtrees_reduced)), jointtree_reduced, dims=1), [[x] for x in keys(jointtree_reduced.lleaves)])
	# tmp = [length(x)!=length(jointtree.lleaves) for x in supra]
	# println("Found $(sum(tmp)) supra-MCCs that are not the root.")
	for (i,sup) in enumerate(supra)
		if length(sup) < threshold_indiv
			# Constructing sub trees
	    	sub_jointtree = node2tree(prunenode(lca([jointtree_reduced.lnodes[x] for x in sup]))[1])
	    	subtrees = Dict{String, Tree}()
	    	for (s,t) in segtrees_reduced
	    	    subtrees[s] = node2tree(prunenode(lca([t.lnodes[x] for x in sup]))[1])
	    	end

	    	# Scores of MCCs
	    	MCCloc = maximal_coherent_clades(collect(values(subtrees))) # Should be the leaves, so I should maybe remove this line
	    	push!(scores_1, prune_mcc_scores(collect(values(subtrees)), sub_jointtree, MCCloc, nmax = threshold_indiv))
	    	# The following ensures `scores_2` has the same length as scores_1
	    	if length(sup) < threshold_pair
				push!(scores_2, prune_mcc_scores_pairs(collect(values(subtrees)), sub_jointtree, MCCloc, nmax = threshold_pair))
			else
				push!(scores_2, Dict{Tuple{String,String},Int64}())
			end
			if length(sup) < threshold_triplets
				push!(scores_3, prune_mcc_scores_triplets(collect(values(subtrees)), sub_jointtree, MCCloc, nmax = threshold_triplets))
			else
				push!(scores_3, Dict{Tuple{String,String,String},Int64}())
			end

		else
			push!(scores_1, Dict{Tuple{String},Int64}())
			push!(scores_2, Dict{Tuple{String,String},Int64}())
			push!(scores_3, Dict{Tuple{String,String,String},Int64}())
		end
	end
	return supra, scores_1, scores_2, scores_3
end

"""
"""
function filter_strains(tofilter_strains, tofilter_mccs ; confidence=[])
	if !isempty(confidence)
		writedlm("StrainsToFilter.txt",cat(tofilter_strains, confidence, dims=2))
	else
		writedlm("StrainsToFilter.txt", tofilter_strains)
	end
	mkpath("outmsas")
	for s in segments
		filter_fasta("../AugurData/results/aligned_h3n2_$(s).fasta", "outmsas/filtered_h3n2_$(s).fasta", tofilter_strains)
	end
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
	_resolve_trees!(segtrees, jointtree; verbose=false)

Delete null branches above internal nodes. Then resolve all trees in `segtrees` using `jointtree`. 
### Notes
Initial null (*i.e.* too short to bear a mutation) branches above internal nodes are not trusted. --> They are set to 0, corresponding internal nodes are removed.  
Trees are resolved using the joint tree, introducing new very small branches. These are trusted because of topological evidence.  
Finally, leaves that have too short of a branch are also stretched to a minimum threshold.  
"""
function _resolve_trees!(segtrees, jointtree; verbose=false)
	verbose && println("\n### RESOLVING ###\n")
	delete_null_branches!(jointtree.root, threshold = 1/length(jointtree.leaves[1].data.sequence)/3)
	jointtree = node2tree(jointtree.root)
	jointtree = remove_internal_singletons(jointtree)
	check_tree(jointtree)
	# Joint tree as a reference
	for (k,t) in segtrees
		verbose && println("$k...")
		L = length(t.leaves[1].data.sequence)
		delete_null_branches!(t.root, threshold = 1/L/3)
		t = node2tree(t.root)
		t = remove_internal_singletons(t)
		t = resolve_trees(t, jointtree, rtau = 1/L/3, verbose=verbose)
		for (k2, t2) in segtrees
			if k!=k2
				t = resolve_trees(t, t2, rtau = 1/L/3, verbose=verbose)
			end
		end
		resolve_null_branches!(t, tau = 1/L/3)
		segtrees[k] = t
	end
end

"""
Write tree, as well as  branch length files for viewing with auspice. 
"""
function write_trees(tdict)
	mkpath("trees")
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
- s1: individual scores for given supra -- s2: pairwise scores for ... Those are dictionaries `s[m] = score`. 
- if `s1[m] == 1` for some MCC `m`, `push!(toremove, m)`. break
"""
function find_suspicious_mccs_irem(supra, mcc_scores_1, mcc_scores_2, mcc_scores_3)
	smcc = Set()
	for (i,sup) in enumerate(supra)
		s1 = mcc_scores_1[i]
		s2 = mcc_scores_2[i]
		s3 = mcc_scores_3[i]
		if !isempty(s1) && findmin([x for x in values(s1)])[1] == 1
			map(x->push!(smcc,x),findall(x->x==1, s1))
		elseif !isempty(s2) && findmin([x for x in values(s2)])[1] == 1
			map(x->(push!(smcc,x[1]);push!(smcc,x[2])),findall(x->x==1, s2))
		elseif !isempty(s3) && findmin([x for x in values(s3)])[1] == 1
			map(x->(push!(smcc,x[1]);push!(smcc,x[2]);push!(smcc,x[3])),findall(x->x==1, s3))
		end
	end
	return smcc
end


end