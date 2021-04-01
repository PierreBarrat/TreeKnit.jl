"""
	mut_crossmap_to_tree!(t::Tree, json::AbstractString, cmkey=:aa_muts)

Add mutation crossmapping `json` file from augur to `t`. Data is added to `n.data.dat[cmkey][s]` for each tree node `n` and segment `s` (read from the json file). 
"""
function mut_crossmap_to_tree!(t::Tree, json::AbstractString, cmkey=:aa_muts)
	cm = JSON3.read(read(json, String), Dict)["nodes"]
	TreeTools.map_dict_to_tree!(t, cm, symbol=true)
	for n in values(t.lnodes)
		for (s,M) in n.data.dat[cmkey]
			n.data.dat[cmkey][s] = [TreeTools.parse_mutation(m) for m in M]
		end
	end
end

"""
	get_mut_dict(t::Tree, s::AbstractString; cmkey="aa_muts")

Construct dictionary mapping mutations of segment `s` in tree `t` to the number of times they appear. For `n` in `t.lnodes`, looks for mutations in `n.data.dat[cmkey][s]`. 
"""
function get_mut_dict(t::Tree, s::AbstractString; cmkey=:aa_muts)
	md = Dict{TreeTools.Mutation,Int64}() 
	for n in values(t.lnodes)
		for m in n.data.dat[cmkey][s]
			md[TreeTools.parse_mutation(m)] = get(md, TreeTools.parse_mutation(m), 0) + 1
		end
	end
	return md
end

"""
	suspicious_branches!(trees::Dict{T,<:Tree}, cmkey=:aa_muts) where T

Mark suspicious branches in `trees` using cross-mapped mutations. Store the number of suspicious mutations for key `s` on the branch above tree node `n` in `n.data.dat[:suspicious_mut][s]`. 
"""
function suspicious_branches!(trees::Dict{T,<:Tree}, cmkey=:aa_muts) where T
	# Set of mutations on branches for segment s on its own tree
	mut_dict = Dict(s=>get_mut_dict(t, s; cmkey=cmkey) for (s,t) in trees)
	# 
	for (s,t) in trees
		for n in values(t.lnodes)
			n.data.dat[:suspicious_mut] = Dict{T, Int64}()	
		end
	end
	# 
	for (sref,tref) in trees, s in Iterators.filter(!=(sref), keys(trees))
		_suspicious_branches!(tref, s, mut_dict, cmkey)
	end
end
function _suspicious_branches!(t::Tree, s, md::Dict, cmkey)
	for n in values(t.lnodes)
		n.data.dat[:suspicious_mut][s] = 0
		for mut in n.data.dat[cmkey][s]
			if !haskey(md[s], mut)
				n.data.dat[:suspicious_mut][s] += 1
			end
		end
	end
	nothing
end

"""
"""
function prune_suspicious_mccs!(trees::Dict, suspmut_key=:suspicious_mut)
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

function _prune_suspicious_mccs!(trees, suspmut_key=:suspicious_mut)
	S = collect(keys(trees))
	# Get mccs
	mccs = maximal_coherent_clades(collect(values(trees)))
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