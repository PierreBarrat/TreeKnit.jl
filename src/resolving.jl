export resolve_trees
export resolve_null_branches!
export is_unresolved_clade

"""
    resolve_trees(t, tref ; label_init=1, rtau = 1. /L/3)

Resolves clades in `t` using clades in `tref`. Newly introduced nodes are assigned a time `rtau`. 
"""
function resolve_trees(t, tref ; rtau = 1. /max(1, length(first(t.lleaves)[2].data.sequence)/3), verbose=false)
    if !share_labels(t, tref)
        error("`resolve_trees` can only be used on trees that share leaf nodes.")
    end       
    tt = deepcopy(t)
    # Getting starting index for node labeling
    label_init = 1
    for n in values(t.nodes)
        if match(r"RESOLVED", n.label)!=nothing && parse(Int64, n.label[10:end]) >= label_init
            label_init = parse(Int64, n.label[10:end]) + 1
        end
    end


    L = length(first(tt.lleaves)[2].data.sequence)
    label_i = label_init
    rcount = 0
    for (i,cl) in enumerate(keys(tref.lleaves))
        if verbose
            print("$i/$(length(tref.lleaves)) -- $cl                                \r")
        end
        # Leaves are always okay
        croot = tref.lleaves[cl]
        clabel = [cl]   
        flag = true
        while flag && !croot.isroot 
            # Go up one clade in `tref`
            nroot = croot.anc
            nlabel = node_leavesclade_labels(nroot) #[x.label for x in node_leavesclade(nroot)]    
            # Three cases
            # i. `nlabel` is a clade in `t` as well.  
            # --> Go up 
            # ii. `nlabel` is an unresolved clade in `t`. 
            # --> Resolve it, and go up
            # iii. `nlabel` is neither an UR clade nor a clade in `t`.
            # --> `flag=false`, go to next leaf. 
            nodes = [tt.lleaves[x] for x in nlabel]
            if isclade(nodes)
            elseif is_unresolved_clade(nodes) # I'm basically calling the same function twice here. Could easily be optimized. 
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
    verbose && println("\n$rcount clades have been resolved.\n")
    return node2tree(tt.root)
end

"""
"""
function resolve_trees(treelist ; label_init=1, verbose=false)
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
        L = max(1, length(first(t.lleaves)[2].data.sequence)
)        # checklist = Dict(k=>false for k in keys(t.lleaves))   
        verbose && println("\n#######################################")
        for (i,cl) in enumerate(keys(t.lleaves))
            verbose && print("$i/$(length(t.lleaves)) -- $cl                                \r")
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
                        elseif is_unresolved_clade(nodes, checkclade=false)
                            # println("Resolving nodes $([nlabel])")
                            make_unresolved_clade!(nodes, label = "RESOLVED_$(label_i)", tau = 1. /L/3)
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
    verbose && println("\n$rcount clades have been resolved.\n")
    return [node2tree(t.root) for t in treelist_]
end

"""
    resolve_null_branches!(t)

Branches shorter than `tau` are set to `tau`. 
"""
function resolve_null_branches!(t::Tree; tau = 1. /max(1, length(first(t.lleaves)[2].data.sequence)/3), internal=false)
    if !internal
        for n in values(t.lleaves)
            if !ismissing(n.data.tau) && n.data.tau < tau
                n.data.tau = tau
            end
        end
    else
        for n in values(t.nodes)
            if !ismissing(n.data.tau) && n.data.tau < tau
                n.data.tau = tau
            end
        end
    end
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
function is_unresolved_clade(node_list ; checkclade=true)
	# (i)
	if checkclade && isclade(node_list)
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