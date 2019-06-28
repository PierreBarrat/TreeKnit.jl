using TreeTools
using FastaIO

function binary_to_nucleotide(infasta::String, outfasta::String)
	fasta = readfasta(infasta)
	mapping = "ACGT-"
	fasta = unique(x->x[1], fasta) # Unique sequence labels
	for (i,(n,s)) in enumerate(fasta)
		fasta[i] = (n,replace(s, ('0','1','2','3','4')=>(x->mapping[parse(Int64,x)+1])))
	end
	writefasta(outfasta, fasta)
end

function add_consensus(infasta::String, outfasta::String ; name="root")
	fasta = readfasta(infasta)
	a = [parse.(Int64, split(x[2],"")) for x in fasta]
	a = round.(Int64, mean(cat(a'..., dims=1),dims=1))
	writefasta(outfasta, [(name, prod("$x" for x in a))])
	writefasta(outfasta, fasta, "a")
end

"""
Assuming we know the *real* label of internal nodes in all trees. 
If two nodes have the same mrca in all trees, they belong to the same MCC.
Only works correctly for two trees right now. 
"""
function realMCCs(trees::Vararg{Tree})
	# Label nodes based on being common in all tree and having a recombination above. 
	# i.e. --> being the root of an MCC!
	if length(trees) > 2 
		@warn "realMCCs only functionnal for 2 trees"
	end
	nd = Dict()
	tref = first(trees)
	for l in keys(tref.lnodes)
		nd[l] = false
		common = prod(haskey(t.lnodes, l) for t in trees)
		if common
			recomb = false
			if |([t.lnodes[l].isroot for t in trees]...)
				nd[l] = true
			elseif length(unique([t.lnodes[l].anc.label for t in trees]))!=1
				alist = [t.lnodes[l].anc for t in trees]
				if length(alist)==2
					a1 = alist[1]
					a2 = alist[2]
					la1 = node_ancestor_list(a1)
					la2 = node_ancestor_list(a2)
					if !in(a1.label, la2) && !in(a2.label, la1)
						nd[l] = true
					end
				else
					@warn "realMCCs only functionnal for 2 trees"
				end
			end
		end
	end

	# Each node n such that nd[n] is true is the root of an MCC. 
	# Go up from each leave unless we encounter one of these nodes
	mcc = Dict(n=>[] for n in findall(nd))
	for l in keys(tref.lleaves)
		cl = l
		while !nd[cl]
			cl = tref.lnodes[cl].anc.label
		end
		push!(mcc[cl], l)
	end
	return mcc
end


"""
"""
function unname_nodes(t::Tree)
	tt = deepcopy(t)
	for n in values(tt.nodes)
		if !n.isleaf
			n.label = ""
		end
	end
	return node2tree(tt.root)
end

