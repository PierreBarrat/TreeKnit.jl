export pmut
export logPoisson
export tree_distancematrix
export tree_distancematrix_null
export match_history
export ref_findrecomb
export tree_match
export _nodes_match!
export likelihoodratio
export sym_likelihoodratio

"""
	pmut(node1::TreeNode, node2::TreeNode, mu, tau)

Probability of mutation from `node1` to `node2` in time `tau` with mutation rate `mu`.  
This is estimated using the Hamming distance between the two sequences as realization of a Poisson distribution with parameter `mu*tau`. 	
"""
function pmut(node1::TreeNode, node2::TreeNode, mu, tau)
    n = hamming(node1.data.sequence, node2.data.sequence)
    return logPoisson(n, tau*mu)
end


"""
	logPoisson(n::Int64, lambda)

Log of the probability of obtaining `n` events in a Poisson distribution of parameter `lambda`.  
Keyval argument `minlambda`: If two segments of two individuals are identical, a tree inference software will conclude that `lambda=0`. Even in a recombinationless case, a mutation might have occured in the other segment. If `lambda=0`, the probability of this mutation is exactly 0, and its significance infinite. 
"""
function logPoisson(n::Int64, lambda; minlambda = 1)
	if lambda > minlambda
    	lp = -lambda + n*log(lambda)
    	for i in 1:n
      	  lp -= log(i)
    	end
    else
    	lp = -1 + n*log(1)
    	for i in 1:n
      	  lp -= log(i)
    	end
    end
    return lp
end

"""
"""
function tree_distancematrix(tree_ref, tree_test, mu)
	nleaves = length(keys(tree_ref.leaves))
	if length(keys(tree_test.leaves)) != nleaves || length(keys(tree_ref.nodes)) != length(keys(tree_test.nodes))
		error("Trees do not have the same number of nodes")
	end
	dmat = zeros(Float64, nleaves, nleaves)
	for i in 1:nleaves
	    dmat[i,i] = 0.
	    for j in (i+1):nleaves
	        tau = node_divtime(tree_ref.leaves[i], tree_ref.leaves[j]) # Time between leaves in the reference tree
	        p1 = pmut(tree_ref.leaves[i], tree_ref.leaves[j], mu, tau) 
	        p2 = pmut(tree_test.leaves[i], tree_test.leaves[j], mu, tau)
	        dmat[i,j] = p2 - p1
	        dmat[j,i] = dmat[i,j]
	    end
	end
	return dmat
end

"""
"""
function tree_distancematrix_null(tree, mu; rep = 2000)
    nleaves = length(keys(tree.leaves))
	pnull = Array{Float64,1}(undef,0)
	for i in 1:nleaves
	    for j in (i+1):nleaves
	        tau = node_divtime(tree.leaves[i], tree.leaves[j]) 
	        n1 = rand(Poisson(mu*tau), rep)
	        n2 = rand(Poisson(mu*tau), rep)
	        append!(pnull, logPoisson.(n1, mu*tau) - logPoisson.(n2,mu*tau))
	    end
	end
	return pnull
end

"""
	likelihoodratio(tree_ref, tree_test, mu)

Compute likelihood ratio value for all pairs of leaves in `tree_ref` and `tree_test`, taking `tree_ref` as a denominator of the ratio.  
Output is a dictionary `Tuple{:label, :label} => value`. 
"""
function likelihoodratio(tree_ref, tree_test, mu)
	out = Dict{Tuple{fieldtype(TreeNode,:label), fieldtype(TreeNode,:label)}, Float64}()
	for (kx,x) in tree_ref.lleaves
		for (ky,y) in tree_ref.lleaves
			if x.label == y.label
				out[(x.label, y.label)] = 0
			else
				tau = node_divtime(x,y)
				p1 = pmut(x,y,mu,tau)
				p2 = pmut(tree_test.lleaves[kx],tree_test.lleaves[ky],mu,tau)
				out[(x.label, y.label)] = p2 - p1
			end
		end
	end
	return out
end


"""
    sym_likelihoodratio(tree1, tree2)

Output is a dictionary `Tuple{:label, :label} => value`. 
"""
function sym_likelihoodratio(tree1::Tree, tree2::Tree)
    out = Dict{Tuple{fieldtype(TreeNode,:label), fieldtype(TreeNode,:label)}, Float64}()
    for (kx,x) in tree1.lleaves
        for (ky,y) in tree1.lleaves
            if x.label == y.label
                out[(x.label, y.label)] = 0
            else
                n1 = hamming(x.data.sequence, y.data.sequence)
                n2 = hamming(tree2.lleaves[kx].data.sequence, tree2.lleaves[ky].data.sequence)
                lml = (n1+n2)/2
                out[(kx,ky)] = sym_likelihoodratio(n1,n2)
            end
        end
    end
    return out
end

function sym_likelihoodratio(n1::Real,n2::Real)
    if n1 == 0 && n2 == 0 
        return 0
    elseif n1 ==0 && n2 != 0
        return -n2*log(2)
    elseif n1 !=0 && n2==0
        return -n1*log(2)
    else
        return (n1+n2)*log((n1+n2)/2) - n1*log(n1) - n2*log(n2)
    end
end



"""
Find set of coherent clades. 
"""
function coherent_clade(tree, rootkey, cmap::Dict{Tuple{String, String}, Bool})
	# List of labels of leavesclade
    cr = map(x->tree.leaves[x].label, tree_leavesclade(tree, rootkey))
    if is_coherent_clade(cr, cmap)
        return [cr]
    end
    cc_list = []
    for c in tree.nodes[rootkey].child
        map(x->push!(cc_list,x), coherent_clade(tree, node_findkey(c,tree), cmap)) # Push all coherent sub-clades into existing coherent clades list 
    end
    return cc_list
end

"""
true if `cladekeys` is a coherent clade, *ie* all pairs verify `cmap[("i","j")] == true`. 
"""
function is_coherent_clade(cladekeys, cmap::Dict{Tuple{String, String}, Bool})
    for li in cladekeys
        for lj in cladekeys
            if !cmap[(li,lj)]
                return false
            end
        end
    end
    return true
end


#####################
##################### SPLITS
#####################

"""
    new_splits(MCCs, t1::Tree, t2::Tree)

Given MCCs and two trees, what are the new splits introduced in either of the trees. 
Return an array `S` of `SplitList` objects, with `S[i][m]` corresponding to new splits in tree `i` from MCC `m`. 

*Note*: Resulting splits are not always unique once constrained to leaves of MCCs! 
"""
function new_splits(MCCs, t1::Tree, t2::Tree)
    splits1 = SplitList(t1)
    splits2 = SplitList(t2)
    # Splits corresponding to each mcc in both trees 
    MCC_splits = splits_in_mccs(MCCs, t1, t2)
    # Are those new or not? 
    _new_splits!(MCC_splits[2], splits1)
    _new_splits!(MCC_splits[1], splits2)
    #
    return MCC_splits[end:-1:1]
end

"""
    new_splits(tref::Tree, MCCs, t::Tree)

What are the splits introduced in `tref` from `t` using consistent clades `MCCs`? 
Return a single `SplitList` object, with unique splits mapped onto the leaves of `tref`. 
"""
function new_splits(tref::Tree, MCCs, t::Tree)
    # Splits in `tref`
    S_ref = SplitList(tref)
    # Splits corresponding to each mcc in tree `t`
    MCC_splits = splits_in_mccs(MCCs, t)
    # Take the new ones only
    _new_splits!(MCC_splits, S_ref)
    # Map them onto leaves of `tref`
    return unique(TreeTools.map_splits_to_tree(MCC_splits, tref), usemask=false)
end

function _new_splits!(MCC_splits, tree_splits)
    eidx = Int64[]
    for (n,S) in enumerate(MCC_splits)
        idx = Int64[]
        for (i,s) in enumerate(S)
            in(s, tree_splits, S.mask) && push!(idx,i)
        end
        deleteat!(S.splits, idx)
        isempty(S.splits) && push!(eidx, n)
    end
    deleteat!(MCC_splits, eidx)
end

"""
Indices of true splits in `S` w.r. to `Sref`. 
"""
function true_splits(S::SplitList, Sref::SplitList, mask=S.mask)
    idx = Int64[]
    for (i,s) in enumerate(S)
        in(s, Sref, mask) && push!(idx, i)
    end
    return idx
end
true_splits(S, tref::Tree, mask=S.mask) = true_splits(S, SplitList(tref), mask)










    


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

