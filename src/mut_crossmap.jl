"""
	function mut_crossmap_to_tree!(t::Tree, json::AbstractString, segment; cmkey=cmmuts selfmutkey = :selfmuts jsonkey=:aa_muts)

Add mutation crossmapping `json` file from augur to `t`. Data is added to `n.data.dat[cmkey][s]` for each tree node `n` and segment `s != segment` (read from the json file). 
Mutations for key `segment` are added to `n.data.dat[selfmutkey]`. 
"""
function mut_crossmap_to_tree!(t::Tree, json::AbstractString, segment; cmkey=:cmmuts, selfmutkey = :selfmuts, jsonkey=:aa_muts, muttype=Char)
	cm = JSON3.read(read(json, String), Dict)["nodes"]
	# Look for `jsonkey` in cm and map it to tree
	TreeTools.map_dict_to_tree!(t, cm, key=jsonkey)
	for n in values(t.lnodes)
		n.data.dat[cmkey] = Dict(s => [TreeTools.parse_mutation(m, muttype) for m in M] for (s,M) in n.data.dat[jsonkey])
		n.data.dat[selfmutkey] = n.data.dat[cmkey][segment]
		delete!(n.data.dat[cmkey], segment)
		delete!(n.data.dat, jsonkey) 
	end
end

"""
	fasta_crossmap_to_trees!(f::Function, trees::Dict, selfkey=:selfseq, cmkey=:cmseq)

Store sequences of several fasta files onto trees. For each pair `(s, t)` of `trees`: 
- Sequences in fasta file `f(s)` are stored on leaves of `t` under `n.data.dat[selfkey]`
- For each other key `σ` of `trees`, sequences in fasta file `f(σ)` are stored on leaves of `t` under `n.data.dat[cmkey][σ]`.
"""
function fasta_crossmap_to_trees!(f::Function, trees::Dict, selfkey=:selfseq, cmkey=:cmseq)
	for (sref,t) in trees
		TreeTools.fasta2tree!(t, f(sref), selfkey)
		for s in Iterators.filter(!=(sref), keys(trees))
			TreeTools.fasta2tree!(t, f(s), (cmkey, s))
		end
	end
end

"""
	function ancestral_sequences!(trees::Dict)
"""
function ancestral_sequences!(trees::Dict; self=true, crossmapped=true)
	for (sref,tref) in trees
		# self sequences and mutations
		self && self_ancestral_states!(tref)
		# crossmapped sequences and mutations
		crossmapped && crossmapped_ancestral_states!(trees, sref)
	end
end

function self_ancestral_states!(tree, key=:selfseq)
	TreeAlgs.Fitch.fitch_mutations!(tree, :selfmuts, key)
	# TreeTools.compute_mutations!(tree, :selfmuts) do n
	# 	n.isleaf ? n.data.dat[key] : n.data.dat[:selfseq_ancestral] 
	# end
end

function crossmapped_ancestral_states!(trees, sref, key=:cmseq)
	tref = trees[sref]
	for s in Iterators.filter(!=(sref), keys(trees))
		TreeAlgs.Fitch.fitch_mutations!(tref, (:cmmuts, s), (key, s))
		# TreeTools.compute_mutations!(tref, (:cmmuts,s)) do n
		# 	n.isleaf ? n.data.dat[key][s] : n.data.dat[:cmseq_ancestral][s] 
		# end
	end
end


"""
	get_mut_dict(f::Function, tree:Tree)

Construct dictionary mapping mutations of segment `s` in tree `t` to the number of times they appear. For `n` in `t.lnodes`, looks for mutations in `f(n)`. 
"""
function get_mut_dict(f::Function, t::Tree)
	md = Dict{TreeTools.Mutation,Int64}() 
	for n in values(t.lnodes)
		for m in f(n)
			md[TreeTools.parse_mutation(m)] = get(md, TreeTools.parse_mutation(m), 0) + 1
		end
	end
	return md
end

"""
    suspicious_mutations(md_ref, md_new ; n = 100)

Following mutations `m` are suspicious: 
- `m` appears in `md_new` but not in `md_ref`
- `m` appears only once in `md_ref`, but multiple times in `md_new`
"""
function suspicious_mutations(md_ref, md_new)
    unique_mut = Set{Tuple{Int64, Int64, Int64}}()
    new_mut = Set{Tuple{Int64, Int64, Int64}}()
    for (k,v) in md_new
        if get(md_ref, k, -1) == -1
            push!(new_mut, k)
        elseif v > 1 && md_ref[k] == 1
            push!(unique_mut, k)
        end
    end
    return (new_mut, unique_mut)
end


"""
	suspicious_branches!(trees::Dict{T,<:Tree}, cmkey=cmmuts where T

Mark suspicious branches in `trees` using cross-mapped mutations. Store the number of suspicious mutations for key `s` on the branch above tree node `n` in `n.data.dat[:suspicious_muts][s]`. 
"""
function suspicious_branches!(trees::Dict{T,<:Tree}, selfmutkey = :selfmuts, cmkey=:cmmuts) where T
	# Set of mutations on branches for segment s on its own tree
	mut_dict = Dict(s=>get_mut_dict(n -> n.data.dat[selfmutkey], t) for (s,t) in trees)
	# 
	for (s,t) in trees
		for n in values(t.lnodes)
			n.data.dat[:suspicious_muts] = Dict{T, Int64}()	
		end
	end
	# 
	for (sref,tref) in trees, s in Iterators.filter(!=(sref), keys(trees))
		_suspicious_branches!(tref, s, mut_dict, cmkey)
	end
end
function _suspicious_branches!(t::Tree, s, md::Dict, cmkey)
	for n in values(t.lnodes)
		n.data.dat[:suspicious_muts][s] = 0

		for mut in n.data.dat[cmkey][s]
			if !haskey(md[s], mut)
				n.data.dat[:suspicious_muts][s] += 1
			end
		end
	end
	nothing
end

"""
"""
function prune_suspicious_mccs!(trees::Dict, suspmut_key=:suspicious_muts)
	if length(trees) > 2
		error("Received $(length(trees)) trees. Only for pairs of trees.")
	end
	pMCCs = []
	pmccs = _prune_suspicious_mccs!(trees, suspmut_key)
	while !isempty(pmccs)
		append!(pMCCs, pmccs)
		pmccs = _prune_suspicious_mccs!(trees, suspmut_key)
	end
	return pMCCs
end

function _prune_suspicious_mccs!(trees, suspmut_key=:suspicious_muts)
	S = collect(keys(trees))
	# Get mccs
	mccs = naive_mccs(collect(values(trees)))
	# For the root of each, see if it has more than 2 suspicious mutations
	idx = Int64[]
	for (i,m) in enumerate(mccs)
		for (sref,t) in trees
			a = lca(t, m)
			for s in Iterators.filter(!=(sref), S) # Just one value normally
				if haskey(a.data.dat,suspmut_key) && a.data.dat[suspmut_key][s] > 1
					@debug "$(a.label) (mcc $i) with $(a.data.dat[suspmut_key][s]) suspicious mutations in $s for tree of $sref"
					@debug "Will prune $m"
					push!(idx, i)
				end
			end
		end
	end
	unique!(idx)
	# Prune those that do
	for i in idx
		for (s,t) in trees
			prunesubtree!(t, mccs[i])
		end
	end
	# Return pruned mccs
	return mccs[idx]
end


