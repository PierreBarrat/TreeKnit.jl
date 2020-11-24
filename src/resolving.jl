export resolve_trees, resolve_trees!, resolve_trees_ref
export resolve_null_branches!
export is_unresolved_clade

"""
    resolve_trees_ref!(t, tref; rtau = :ref verbose=false)

Resolve clades in place in `t` using clades in `tref`. Newly introduced nodes are assigned a time `rtau`. If `rtau==:ref`, the time found in `tref` is used. 
"""
function resolve_trees_ref!(t, tref; rtau = :ref, verbose=false)
    rcount = _resolve_trees_ref!(t, tref, rtau, verbose)
    verbose && println("$rcount clades have been resolved.")
    return nothing
end
"""
    resolve_trees_ref(t, tref; rtau = :ref verbose=false)

Resolve clades in `t` using clades in `tref`. Return a new tree, `t` is unchanged. Newly introduced nodes are assigned a time `rtau`. If `rtau==:ref`, the time found in `tref` is used. 
"""
function resolve_trees_ref(t, tref; rtau = :ref, verbose=false)
    tc = deepcopy(t)
    resolve_trees_ref!(tc, tref, rtau=rtau, verbose=verbose)
    return tc
end

"""
    resolve_trees(treelist::Vararg{Tree}; rtau=:ref, verbose=false)
    resolve_trees!(treelist::Vararg{Tree}; rtau=:ref, verbose=false)

Use each tree of treelist in turn as a reference. While resolutions are found keep cycling. 
"""
function resolve_trees!(treelist::Vararg{Tree}; rtau=:ref, verbose=false)
    found_resolution = true
    while found_resolution
        flag = false
        for (iref, tref) in enumerate(treelist)
            for (i,t) in Iterators.filter(x->x[1]!=iref, enumerate(treelist))
                verbose && println("\nRef: $iref - object: $i")
                rcount = _resolve_trees_ref!(t, tref, rtau)
                rcount > 0 && (flag = true)
                verbose && rcount > 0 && println("Resolved $rcount clades in tree $i with tree $iref as reference.")
            end
        end
        found_resolution = flag
    end
    return treelist
end
function resolve_trees(treelist::Vararg{Tree}; rtau=:ref, verbose=false)
    ntl = deepcopy(treelist)
    resolve_trees!(ntl..., rtau=rtau, verbose=verbose)
    return ntl
end

function _resolve_trees_ref!(t, tref, rtau)
    
    !share_labels(t, tref) && @error "Can only resolve trees that share leaf nodes."
    
    # Get starting index for node labeling
    label_init = 1
    for n in values(t.lnodes)
        if match(r"RESOLVED", n.label)!=nothing && parse(Int64, n.label[10:end]) >= label_init
            label_init = parse(Int64, n.label[10:end]) + 1
        end
    end
    #
    label_i = label_init
    rcount = 0
    for (i,cl) in enumerate(keys(tref.lleaves))
        nroot = tref.lleaves[cl]
        while !nroot.isroot 
            # Go up one clade in `tref`
            nroot = nroot.anc
            nlabel = node_leavesclade_labels(nroot) # Clade in `tref`
            # Three cases
            # i. `nlabel` is a clade in `t` as well.  
            # --> Go up 
            # ii. `nlabel` is an unresolved clade in `t`. 
            # --> Resolve it, and go up
            # iii. Neither
            # --> Break, go to next leaf. 
            nodes = [t.lleaves[x] for x in nlabel]
            if !isclade(nodes) && is_unresolved_clade(nodes) 
                tau = (rtau == :ref) ? nroot.data.tau : rtau
                make_unresolved_clade!(nodes, label = "RESOLVED_$(label_i)", tau = tau, safe=false)
                rcount += 1
                label_i += 1
            else
                break
            end
        end
    end
    node2tree!(t, t.root)
    return rcount
end

# function resolve_trees(treelist::Vararg{Tree})

# """
# """
# function resolve_trees!(treelist; label_init=1, verbose=false, τ=missing)

#     sh = mapreduce(t->share_labels(t[1],t[2]), *, zip(treelist[1:end-1], treelist[2:end]))
#     !sh && @error "Can only resolve trees that share leaf nodes."

#     for t in treelist
#         for n in values(t.lnodes)
#             if match(r"RESOLVED", n.label)!=nothing && parse(Int64, n.label[10:end]) >= label_init
#                 label_init = parse(Int64, n.label[10:end]) + 1
#             end
#         end
#     end

#     label_i = label_init
#     rcount = 0
#     # Use clades of each tree `t` to resolve other trees
#     for t in treelist    
#         L = max(1, length(first(t.lleaves)[2].data.sequence))
#         ismissing(τ) ? τe = 1/L/3 : τe = τ
#         # checklist = Dict(k=>false for k in keys(t.lleaves))   
#         # verbose && println("\n#######################################")
#         for (i,cl) in enumerate(keys(t.lleaves))
#             # verbose && print("$i/$(length(t.lleaves)) -- $cl                                \r")
#             # println("\nStarting from $(cl)")
#             croot = t.lleaves[cl]
#             clabel = [cl]
#             flag = true
#             while flag && !croot.isroot 
#                 # Go up one clade in current tree `t`
#                 # `nlabel` are the labels of nodes in that clade
#                 nroot = croot.anc
#                 nlabel = [x.label for x in node_leavesclade(nroot)]
#                 # Three cases
#                 # i. `nlabel` is a clade in all trees. 
#                 # --> Go up in all trees 
#                 # ii. `nlabel` is a clade in some trees, an unresolved one some of the other trees
#                 # --> Resolve it in the other trees, and go up
#                 # iii. `nlabel` is neither an UR clade nor a clade in some trees.
#                 # --> `flag=false`, go to next leaf. 
#                 # println("Checking clade based at $(nroot.label).")
#                 # println("Labels are $(nlabel)")
#                 for ct in treelist
#                         nodes = [ct.lleaves[x] for x in nlabel]
#                         if isclade(nodes)
#                         elseif is_unresolved_clade(nodes, checkclade=false)
#                             # println("Resolving nodes $([nlabel])")
#                             make_unresolved_clade!(nodes, label = "RESOLVED_$(label_i)", tau = τe)
#                             rcount += 1
#                             label_i += 1
#                         else
#                             flag = false
#                         end
#                 end
#                 if flag
#                     croot = nroot
#                     clabel = nlabel
#                 end
#             end
#         end
#     end
#     verbose && println("\n$rcount clades have been resolved.\n")
#     for (i,t) in enumerate(treelist)
#         node2tree!(t, t.root)
#     end
#     return nothing
# end


# """
# """
# function resolve_trees(treelist ; label_init=1, verbose=false, τ=missing)
#     flag = true
#     for t1 in treelist
#         for t2 in treelist
#             flag *= share_labels(t1,t2)
#         end
#     end
#     for t in treelist
#         for n in values(t.lnodes)
#             if match(r"RESOLVED", n.label)!=nothing && parse(Int64, n.label[10:end]) >= label_init
#                 label_init = parse(Int64, n.label[10:end]) + 1
#             end
#         end
#     end
#     if !flag
#         error("`resolve_trees` can only be used on trees that share leaf nodes.")
#     end

#     treelist_ = deepcopy(treelist)
#     label_i = label_init
#     rcount = 0
#     # Use clades of each tree `t` to resolve other trees
#     for t in treelist_    
#         L = max(1, length(first(t.lleaves)[2].data.sequence))
#         ismissing(τ) ? τe = 1/L/3 : τe = τ
#         # checklist = Dict(k=>false for k in keys(t.lleaves))   
#         # verbose && println("\n#######################################")
#         for (i,cl) in enumerate(keys(t.lleaves))
#             # verbose && print("$i/$(length(t.lleaves)) -- $cl                                \r")
#             # println("\nStarting from $(cl)")
#             croot = t.lleaves[cl]
#             clabel = [cl]
#             flag = true
#             while flag && !croot.isroot 
#                 # Go up one clade in current tree `t`
#                 # `nlabel` are the labels of nodes in that clade
#                 nroot = croot.anc
#                 nlabel = [x.label for x in node_leavesclade(nroot)]
#                 # Three cases
#                 # i. `nlabel` is a clade in all trees. 
#                 # --> Go up in all trees 
#                 # ii. `nlabel` is a clade in some trees, an unresolved one some of the other trees
#                 # --> Resolve it in the other trees, and go up
#                 # iii. `nlabel` is neither an UR clade nor a clade in some trees.
#                 # --> `flag=false`, go to next leaf. 
#                 # println("Checking clade based at $(nroot.label).")
#                 # println("Labels are $(nlabel)")
#                 for ct in treelist_
#                         nodes = [ct.lleaves[x] for x in nlabel]
#                         if isclade(nodes)
#                         elseif is_unresolved_clade(nodes, checkclade=false)
#                             # println("Resolving nodes $([nlabel])")
#                             make_unresolved_clade!(nodes, label = "RESOLVED_$(label_i)", tau = τe)
#                             rcount += 1
#                             label_i += 1
#                         else
#                             flag = false
#                         end
#                 end
#                 if flag
#                     croot = nroot
#                     clabel = nlabel
#                 end
#             end
#         end
#     end
#     verbose && println("\n$rcount clades have been resolved.\n")
#     return [node2tree(t.root) for t in treelist_]
# end



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
(i) If `node_list` forms a clade, return `false`. Else:  

(ii) Let `MRCA` be the common ancestor of all members of `node_list`. For each `n` in `node_list`, find the ancestor just below `MRCA`, *i.e.* the ancestor just before all members of `node_list` join. The set of these ancestors is `A`.  

(iii) For each member `a` of `A`, construct the corresponding clade. All nodes in this clade should also belong to `node_list`. Why? All members of `A` can be grouped into a clade, and this clade would be exactly `node_list`. 

OLD(iii) A copy of the tree is created. Each element in `A` is pruned of the copy and grafted onto an empty root `r`. If the leaves of this new tree correspond exactly to `node_list`, then `node_list` is an unresolved clade, and `r` is the root of this clade. 

"""
function is_unresolved_clade(node_list ; checkclade=true)
	# (i)
	if checkclade && isclade(node_list)
		return false
	end
    
    label_list = [x.label for x in node_list]
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
        # if !issubset(Set(x.label for x in POTleaves(a)), label_list)
        for n in POTleaves(a)
            if !in(n.label, label_list)
                return false
            end
        end
	end
    return true
end


"""
    make_unresolved_clade!(nodelist)

Resolve the unresolved clade defined by `nodelist` by creating a new internal node `r` in the tree. Return `r`, 
"""
function make_unresolved_clade!(node_list ; label="", tau = 0., safe=true)

    if safe && !is_unresolved_clade(node_list)
        error("Input is not an unresolved clade")
    end


    mrca = lca(node_list)
    checklist = Dict((x.label)=>false for x in node_list)
    r = TreeNode()
    for n in node_list
        if !checklist[n.label]
            a = n
            while (a.anc != mrca) # && !(a.isroot)
                a = a.anc
            end
            map(x->checklist[x.label]=true, POTleaves(a))
            prunenode!(a)
            graftnode!(r,a)
        end
    end

    r.label = label
    graftnode!(mrca, r, tau = tau)    
    return r
end