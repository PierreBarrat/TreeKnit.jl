"""
1. Resolve trees
2. Find MCCs
3. Infer common subtrees for MCCs using joint sequences
"""
function prunetrees_mut(_segtrees, jointmsa; verbose = true, cwd="$(pwd())")
	segtrees = deepcopy(_segtrees)

	## Resolving, finding MCCs, and adjusting branch length to that of the joint tree
	# resolving
	verbose && println("Resolving trees...")
	_resolve_trees!(segtrees);
	# MCCs
	MCC = maximal_coherent_clades(collect(values(segtrees)))
	verbose && println("Found $(length(MCC)) MCCs of average size $(mean([length(x) for x in MCC]))")
	# Adjusting branch length in MCCs
	mcc_names = name_mcc_clades!(collect(values(segtrees)), MCC)
	verbose && println("Adjusting branch length in MCCs: inferring $(length(MCC)) subtrees using full genome...")
	_adjust_branchlength!(segtrees, jointmsa, MCC)

	## Computing ancestral states
	verbose && println("Cross-mapping sequences on trees...")
	crossmapping(segtrees, cwd)

	## Counting mutations
	map(s->fasta2tree!(segtrees[s], "$cwd/CrossMapping/ancestral_tree_$(s)_aln_$(s)/ancestral_sequences.fasta"), segments)
	segmd_ref = Dict(s=>make_mutdict!(segtrees[s]) for s in segments)
	crossmuts = compute_crossmuts(segtrees, MCC, segmd_ref, cwd)

	## Finding suspicious MCCs
	fmcc = find_suspicious_mccs(crossmuts, segtrees)
	# iroot = findall(x->x==root_strain, tofilter)
	for (i,m) in enumerate(fmcc)
		if in(root_strain, mcc_names[m])
			println(i)
			splice!(fmcc,i)
		end
	end
	# return tofilter, fmcc, confidence
	return segtrees, fmcc, mcc_names
end



"""
	_adjust_branchlength(segtrees, jointmsa, MCC)
"""
function _adjust_branchlength!(segtrees, jointmsa, MCC)
	for (s,t) in segtrees
		segtrees[s] = _adjust_branchlength(t, jointmsa, MCC)
	end
end

"""
	_adjust_branchlength(tree::Tree, jointmsa, MCC)

Copy `tree`. For all `m` in MCC: 
1. Prune MRCA of `m` from `tree`, giving `sub`
2. Scale branches of `sub` by calling Interfacing.scalebranches
3. Regraft `sub` at the right position
"""
function _adjust_branchlength(tree::Tree, jointmsa, MCC)
	t = deepcopy(tree)
	for m in MCC
		if length(m) > 2
			# Subtree to resclae
			sub, r = prunesubtree!(t, m);
			tau = sub.root.data.tau
			rootlabel = sub.root.label
			# Rescaling
			sub = RecombTools.Interfacing.scalebranches(sub, jointmsa);
			# Rerooting 
			root = node_findlabel(rootlabel, sub.root)
			TreeTools.reroot!(root)
			root.data.tau = tau
			# Regrafting
			graftnode!(r, root)
			t = node2tree(t.root);
		end
	end
	return t
end


"""
"""
function crossmapping(segtrees, cwd::String)
	mkpath("$cwd/CrossMapping")
	for (s,t) in segtrees
		write_newick("$cwd/CrossMapping/tree_$(s)_adjusted.nwk", t)
	end
	map(k->write_fasta("$cwd/CrossMapping/aligned_simple_h3n2_$(k).fasta", segtrees[k], internal=false), collect(keys(segtrees)))

	for s1 in segments
		for s2 in segments
			RecombTools.Interfacing.run_treetime(verbose=false, aln="$cwd/CrossMapping/aligned_simple_h3n2_$(s1).fasta", tree="$cwd/CrossMapping/tree_$(s2)_adjusted.nwk", out="$cwd/CrossMapping/ancestral_tree_$(s2)_aln_$(s1)")
		end
	end
end

"""
Construct a dictionary of the form `"mcc label"=>CrossMutations`. 
Expects the directory CrossMapping to exist, containing ancestral sequences for all combinations of segments and trees. 
"""
function compute_crossmuts(segtrees, MCC, segmd_ref, cwd)
	# Building a [tree, aln] dictionary of mutations
	crossmuts_all = Dict()
	for a in keys(segtrees)
	    for t in keys(segtrees)
	        crossmuts_all[t,a] = parse_nexus("$cwd/CrossMapping/ancestral_tree_$(t)_aln_$(a)/annotated_tree.nexus")
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
	tofilter = Array{Any,2}(undef, 0, 2)
	confidence = []
	for (l,cm) in crossmuts # loop on mccs
	    if sum(values(cm.suspicious))>0
	        # push!(fmcc, node_leavesclade_labels(segtrees[segments[1]].lnodes[l]))
	        push!(fmcc, l)
	        # push!(tofilter, [x.label for x in node_leavesclade(segtrees[segments[1]].lnodes[l])]...)
	        # strains = node_leavesclade_labels(segtrees[segments[1]].lnodes[l])
	        # tofilter = [tofilter ; cat(strains, repeat([mcc_idx()], length(strains)), dims=2)]
        	# push!(confidence, repeat([sum(collect(values(cm.suspicious)))], length(node_leavesclade(segtrees[segments[1]].lnodes[l])))...)
	    end
	end
	return fmcc
end