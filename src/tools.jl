export pmut
export logPoisson
export tree_distancematrix
export tree_distancematrix_null
export coherent_clade
export is_coherent_clade
export match_history
export ref_findrecomb
export tree_match
export _nodes_match!
export clade_maptonodes
export likelihoodratio
export sym_likelihoodratio
export is_common_clade
export maximal_coherent_clades, maximal_coherent_clades_old
export is_unresolved_clade, is_coherent_clade_nodelist
export is_coherent_clade
export make_unresolved_clade!
export resolve_trees
export resolve_trees_new
export name_mcc_clades!
export suspicious_mutations
export resolve_null_branches!

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

function sym_likelihoodratio(n1::Int64,n2::Int64)
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

"""
    resolve_trees(t, tref ; label_init=1, rtau = 1. /L/4)

Resolves clades in `t` using clades in `tref`. Newly introduced nodes are assigned a time `rtau`. 
"""
function resolve_trees(t, tref ; label_init=1, rtau = 1. /length(t.leaves[1].data.sequence)/4)
    if !share_labels(t, tref)
        error("`resolve_trees` can only be used on trees that share leaf nodes.")
    end       
    tt = deepcopy(t)
    L = length(tt.leaves[1].data.sequence)
    label_i = label_init
    rcount = 0
    for (i,cl) in enumerate(keys(tref.lleaves))
        print("$i/$(length(tref.lleaves)) -- $cl                                \r")
        # Leaves are always okay
        croot = tref.lleaves[cl]
        clabel = [cl]   
        flag = true
        while flag && !croot.isroot 
            # Go up one clade in `tref`
            nroot = croot.anc
            nlabel = [x.label for x in node_leavesclade(nroot)]    
            # Three cases
            # i. `nlabel` is a clade in `t` as well.  
            # --> Go up 
            # ii. `nlabel` is an unresolved clade in `t`. 
            # --> Resolve it, and go up
            # iii. `nlabel` is neither an UR clade nor a clade in `t`.
            # --> `flag=false`, go to next leaf. 
            nodes = [tt.lleaves[x] for x in nlabel]
            if isclade(nodes)
            elseif is_unresolved_clade(nodes)
                make_unresolved_clade!(nodes, label = "RESOLVED_$(label_i)", tau = rtau)
                rcount += 1
                label_i += 1
            else
                flag = false
            end
            if flag
                croot = nroot
                clabel = nlabel
            end
        end
    end
    println("\n$rcount clades have been resolved.\n")
    return node2tree(tt.root)
end

"""
"""
function resolve_trees(treelist ; label_init=1)
    flag = true
    for t1 in treelist
        for t2 in treelist
            flag *= share_labels(t1,t2)
        end
    end
    if !flag
        error("`resolve_trees` can only be used on trees that share leaf nodes.")
    end

    treelist_ = deepcopy(treelist)
    label_i = label_init
    rcount = 0
    # Use clades of each tree `t` to resolve other trees
    for t in treelist_    
        L = length(t.leaves[1].data.sequence)
        # checklist = Dict(k=>false for k in keys(t.lleaves))   
        println("\n#######################################")
        for (i,cl) in enumerate(keys(t.lleaves))
            print("$i/$(length(t.lleaves)) -- $cl                                \r")
            # println("\nStarting from $(cl)")
            croot = t.lleaves[cl]
            clabel = [cl]
            flag = true
            while flag && !croot.isroot 
                # Go up one clade in current tree `t`
                # `nlabel` are the labels of nodes in that clade
                nroot = croot.anc
                nlabel = [x.label for x in node_leavesclade(nroot)]
                # Three cases
                # i. `nlabel` is a clade in all trees. 
                # --> Go up in all trees 
                # ii. `nlabel` is a clade in some trees, an unresolved one some of the other trees
                # --> Resolve it in the other trees, and go up
                # iii. `nlabel` is neither an UR clade nor a clade in some trees.
                # --> `flag=false`, go to next leaf. 
                # println("Checking clade based at $(nroot.label).")
                # println("Labels are $(nlabel)")
                for ct in treelist_
                        nodes = [ct.lleaves[x] for x in nlabel]
                        if isclade(nodes)
                        elseif is_unresolved_clade(nodes)
                            # println("Resolving nodes $([nlabel])")
                            make_unresolved_clade!(nodes, label = "RESOLVED_$(label_i)", tau = 1. /L/4)
                            rcount += 1
                            label_i += 1
                        else
                            flag = false
                        end
                end
                if flag
                    croot = nroot
                    clabel = nlabel
                end
            end
        end
    end
    println("\n$rcount clades have been resolved.\n")
    return [node2tree(t.root) for t in treelist_]
end

"""
    resolve_null_branches!(t)
"""
function resolve_null_branches!(t::Tree; tau = 1. /length(t.leaves[1].data.sequence)/4, internal=false)
    if !internal
        for n in values(t.leaves)
            if n.data.tau < tau
                n.data.tau = tau
            end
        end
    else
        for n in values(t.nodes)
            if n.data.tau < tau
                n.data.tau = tau
            end
        end
    end
end


"""
    maximal_coherent_clades(treelist)

Find sets of nodes which are: 
- clades in all trees of `treelist`,
- all subclades of nodes are clades in all trees of `treelist` (both of these properties define consistency),
- maximal: adding a node to a set results it in not being a clade in at least one of the trees. 
All the trees of `treelist` should share the same leaf nodes.  

# Note
In this version, the function does not attempt to
- Resolve clades. Since we should already be resolving clade using the information of all segments, resolving them here just makes the code more complex
- Increasing MCCs by adding children of multiforcations one by one. I wish to keep this code as basic as possible: it should just find regions of perfect topologic compatibility in all trees of `treelist`. The rest can use another function, maybe a `combine_mcc` one. 
"""
function maximal_coherent_clades(treelist)
    # Checking that trees have the same label for leaf nodes
    flag = true
    for t1 in treelist
        for t2 in treelist
            flag *= share_labels(t1,t2)
        end
    end
    if !flag
        error("`maximal_common_clades` can only be used on trees that share leaf nodes.")
    end
    # List of already visited nodes
    treelist_ = deepcopy(treelist)
    t = treelist_[1]
    checklist = Dict(k=>false for k in keys(t.lleaves))
    # Explore one leave at a time
    mc_clades = []
    global i = 1;
    println("Hello")
    for (cl,v) in checklist
        print("$i/$(length(t.lleaves)) -- $cl  --     Found $(sum(values(checklist)))/$(length(t.lleaves))            \r")
        if !v # If leave not already visited
            # We're going to go up in all trees at the same time
            global croot = map(x->x.lleaves[cl], treelist_) 
            global clabel = [cl]
            # Initial individual, always a common clade in all trees since it's a leaf. 
            flag = true
            while flag && prod([!x.isroot for x in croot])
                nroot = [x.anc for x in croot] # Ancestors of current label in all trees
                # Each element of `nroot` defines a set of labels. There are two possibilites 
                # (i) Those sets of labels match. In this case, we have a potential consistent clade. To check further, call `is_coherent_clade`. 
                # (ii) Otherwise, the topology of trees in `treelist` is inconsistent above `croot`. `croot` is an MCC, break.  
                nlabel = [Set(x.label for x in node_leavesclade(r)) for r in nroot] # List of sets of labels
                if prod([nlabel[i]==nlabel[1] for i in 1:length(nlabel)])
                    nclade = [[node_findlabel(l, r) for l in nlabel[1]] for r in nroot]
                    if prod([is_coherent_clade_nodelist(c, treelist_) for c in nclade])
                        croot = [lca(c) for c in nclade]
                        clabel = nlabel[1]
                    else
                        flag = false
                    end
                else
                    flag = false
                end
            end

            ###
            # clabel = map(x->x.label, cclade)
            map(x->checklist[x]=true, [c for c in clabel])
            push!(mc_clades, [c for c in clabel])
        end
        global i+=1
    end
    return mc_clades
end


"""
    maximal_coherent_clades_old(treelist)

Find sets of nodes which are: 
- clades in all trees of `treelist`,
- all subclades of nodes are clades in all trees of `treelist` (both of these properties define consistency),
- maximal: adding a node to a set results it in not being a clade in at least one of the trees. 
All the trees of `treelist` should share the same leaf nodes.  
"""
function maximal_coherent_clades_old(treelist)
    # Checking that trees have the same label for leaf nodes
    flag = true
    for t1 in treelist
        for t2 in treelist
            flag *= share_labels(t1,t2)
        end
    end
    if !flag
        error("`maximal_common_clades` can only be used on trees that share leaf nodes.")
    end
    # List of already visited nodes
    treelist_ = deepcopy(treelist)
    t = treelist_[1]
    checklist = Dict(k=>false for k in keys(t.lleaves))
    # Explore one leave at a time
    mc_clades = []
    global i = 1;
    for (cl,v) in checklist
        print("$i/$(length(t.lleaves)) -- $cl  --     Found $(sum(values(checklist)))/$(length(t.lleaves))            \r")
        if !v # If leave not already visited
            ###
            # We're going to go up in all trees at the same time
            global croot = map(x->x.lleaves[cl], treelist_) 
            global clabel = [cl]
            # Initial individual, always a common clade in all trees since it's a leaf. 
            flag = true
            while flag && prod([!x.isroot for x in croot])
                nroot = [x.anc for x in croot] # Ancestors of current label in all trees
                # Each element of `nroot` defines a set of labels. There are four possibilites
                # (i) The intersection of those labels is equal to the labels defined by `clabel` (*i.e.* the clade defined by `croot`). This means that the trees completely diverge above `croot` because of recombination. The current clade is an MCC, break.  
                # (ii) The intersection of those labels defines a coherent resolved clade in *all* trees. In this case, assign the root of this clade to croot, for all trees, and iterate
                # (iii) The intersection of those labels defines a resolved clade in part of the trees, and an unresolved clade in others. In this case, introduce an internal node to resolve the clades in the non resolved trees. Assign the roots of the already existing and newly formed clades to croot. 
                # (iv) None of the above: this means  there has been some recombination. We must now try to build a maximally UR clade, pruning the recombinant parts. 
                # Try to augment the current maximal clade defined by croot by attempting to add each child of each element of `nroot` (one by one). When no more child can be added, the current clade is an MCC, break.  
                nlabel = intersect(map(r->[x.label for x in node_leavesclade(r)],nroot)...)  
                nclade = [[node_findlabel(l, r) for l in nlabel] for r in nroot]
                if nlabel == clabel
                    flag = false
                elseif prod([is_coherent_clade_nodelist(c, treelist_) for c in nclade])
                    croot = [lca(c) for c in nclade]
                    clabel = nlabel
                else
                    flag = false
                    for r in nroot
                        for c in r.child
                            tlabel = intersect(union(clabel, [x.label for x in node_leavesclade(c)]), nlabel)
                            tclade = [node_findlabel(l, treelist_[1].root) for l in tlabel]
                            if tlabel != clabel && is_coherent_clade_nodelist(tclade, treelist_)
                                clabel = tlabel
                            end
                        end
                    end
                end
            end

            ###
            # clabel = map(x->x.label, cclade)
            map(x->checklist[x]=true, clabel)
            push!(mc_clades, clabel)
        end
        global i+=1
    end
    return mc_clades
end


            

"""
Check whether a coherent clade can be built out of `nodelist`. 
i. Check that `nodelist` is either a clade or an unresolved clade. 
ii. Find the common ancestor to `nodelist`. Create it if it is unresolved. 
iii. Check whether this common ancestor is a coherent clade  
All members of `nodelist` should be leaves.
"""
function is_coherent_clade_nodelist(nodelist::Array{TreeNode,1}, treelist)
    if !mapreduce(x->x.isleaf, *, nodelist)
        error("All nodes in `nodelist` should be leaves.")
    end
    if length(nodelist)==1
        return is_coherent_clade(nodelist[1], treelist)
    end
    
    if isclade(nodelist)
        A = lca(nodelist)
    # elseif ((f, A) = is_unresolved_clade(nodelist))[1] ## TO RECODE
    else
        return false
    end
    return is_coherent_clade(A, treelist)
end

"""
Check whether the clade defined by `node` is a coherent clade: 
- it is a clade in all trees of `treelist`
- it is a clade at all levels, *i.e.* for all children `c` of `node`, `is_coherent_clade(c, treelist)` is true
"""
function is_coherent_clade(node::TreeNode, treelist)
    # If it's a leaf, it's a coherent clade
    if node.isleaf
        return true
    end
    # Is `node` a common clade to all trees? 
    cl = map(x->x.label, node_leavesclade(node))
    if !is_common_clade(cl, treelist)
        return false
    end
    # If yes, are all of its children coherent clades? 
    for c in node.child
        if !is_coherent_clade(c, treelist)
            return false
        end
    end
    return true
end

"""
"""
function is_common_clade(label_list, treelist)
    out = true
    for tree in treelist
        nodelist = map(x->tree.lleaves[x], label_list)
        out *= isclade(nodelist)
    end
    return out
end


"""
"""
function clade_maptonodes(cmap)
	cl = []
	l = 1
	c = 1
	while l!=0
		tmp = findall(x->x==c, cmap)
		l = length(tmp)
		if l >0
			push!(cl, tmp)
		end
		c+=1
	end
	return cl
end

"""
	is_unresolved_clade(node_list)

`node_list` is an unresolved clade on one of two conditions:
- `node_list` defines a clade.
- `node_list` could define a clade by creation of only one branch in the tree, adding an exclusive common ancestor to all of its members. 

# Algorithm
(i) If `node_list` forms a clade, return `true`. Else:  

(ii) Let `MRCA` be the common ancestor of all members of `node_list`. For each `n` in `node_list`, find the ancestor just below `MRCA`, *i.e.* the ancestor just before all members of `node_list` join. The set of these ancestors is `A`.  

(iii) For each member `a` of `A`, construct the corresponding clade. All nodes in this clade should also belong to `node_list`. Why? Because this means we could group all members of `A` into a clade and recover exactly `node_list`. 

OLD(iii) A copy of the tree is created. Each element in `A` is pruned of the copy and grafted onto an empty root `r`. If the leaves of this new tree correspond exactly to `node_list`, then `node_list` is an unresolved clade, and `r` is the root of this clade. 

"""
function is_unresolved_clade(node_list)
	# (i)
	if isclade(node_list)
		return true
	end
    
    label_list = Set(x.label for x in node_list)
    mrca = lca(node_list)
    for (i,n) in enumerate(node_list)
	   # (ii)
		a = n
		while a.anc != mrca
            if a.isroot
                println(n.label)
            end
			a = a.anc
		end
        # (iii)
        if !issubset(Set(x.label for x in node_leavesclade(a)), label_list)
            return false
        end
	end
    return true

    # OLD
	# r = TreeNode()
	# for a in A
	# 	aa = prunenode(a)
	# 	graftnode!(r, aa)	
	# end
	# r_clade = Set(x.label for x in node_leavesclade(r))
	# return (Set(x.label for x in node_list) == r_clade), r
end


"""
"""
function name_mcc_clades!(treelist, MCC ; label_init = 1)
    roots = [t.root for t in treelist]
    tt = label_init
    for (i, m) in enumerate(MCC)
        ft = false
        for t in roots
            # Find said clade
            tclade = [node_findlabel(l,t) for l in m] 
            # Introduce new internal node
            if isclade(tclade)
                r = lca(tclade)
                if r.isleaf
                    r.label = "MCC_$(r.label)"
                else
                    r.label = "MCC_$tt"
                    ft = true
                end
            else
                @warn "Should be either a clade or an unresolved clade."
            end
        end
        if ft
            tt += 1
        end
    end
    for t in treelist
        tt = node2tree(t.root)
        t.lleaves = tt.lleaves
        t.lnodes = tt.lnodes
    end
end

    
"""
    make_unresolved_clade!(nodelist)

Resolve the unresolved clade defined by `nodelist` by creating a new internal node `r` in the tree. Return `r`, 
"""
function make_unresolved_clade!(node_list ; label="", tau = 0.)
    if isclade(node_list)
        return lca(node_list)
    end
    if !is_unresolved_clade(node_list)
        error("`node_list` is not an unresolved clade")
    end


    mrca = lca(node_list)
    checklist = Dict((x.label)=>false for x in node_list)
    r = TreeNode()
    for n in node_list
        if !checklist[n.label]
            a = n
            while (a.anc != mrca) && !(a.isroot)
                a = a.anc
            end
            map(x->checklist[x.label]=true, node_leavesclade(a))
            prunenode!(a)
            graftnode!(r,a)
        end
    end

    r.label = label
    graftnode!(mrca, r, tau = tau)    
    return r

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

