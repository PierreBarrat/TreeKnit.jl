export resolve_trees, resolve_trees!, resolve_trees_ref
export resolve_null_branches!
export is_unresolved_clade

"""
    resolve_trees_ref!(t, tref; rtau = :ref, verbose=false)

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

function resolve_trees_joint!(treelist::Vararg{Tree}; verbose=false)
    found_resolution = true
    while found_resolution
        flag = false
        for (itar, targettree) in enumerate(treelist)
            tref = [x[2] for x in Iterators.filter(x->x[1]!=itar, enumerate(treelist))]
            verbose && println("\nTarget: $itar")
            rcount = resolve_trees_target!(targettree, tref..., verbose=verbose)
            rcount > 0 && (flag = true)
            verbose && rcount > 0 && println("Resolved $rcount clades in tree $itar.")
        end
        found_resolution = flag
    end
end
"""
    resolve_trees_target!(t::Tree, tref::Vararg{Tree}; verbose=false)

Resolve target tree `t` using trees in `tref`. Only non-conflicting resolutions are kept. Time of introduced branches are the average of the times of roots of clades in guiding trees. 
"""
function resolve_trees_target!(t::Tree, tref::Vararg{Tree}; verbose=false)
    # Get starting index for node labeling
    label_i = 1
    for n in values(t.lnodes)
        if match(r"RESOLVED", n.label)!=nothing && parse(Int64, n.label[10:end]) >= label_i
            label_i = parse(Int64, n.label[10:end]) + 1
        end
    end
    # 
    resolvable_clades = _resolve_trees_target(t, tref..., verbose=verbose)
    tau = 0.
    Z = 0
    for c in resolvable_clades
        # Getting an appropriate time
        for tr in tref
            clade = [tr.lleaves[x] for x in c]
            if isclade(clade)
                a = lca(clade)
                !ismissing(a) && (tau+=a.data.tau; Z += 1)
            end
        end
        Z == 0. ? (tau = missing) : (tau/=Z)
        # Resolving
        make_unresolved_clade!([t.lleaves[x] for x in c], label = "RESOLVED_$(label_i)", tau=tau, safe=false)
        label_i += 1
    end
    node2tree!(t, t.root)
    check_tree(t)
    return length(resolvable_clades)
end
"""
    _resolve_trees_target(t::Tree, tref::Vararg{Tree}; verbose=false)

Return resolvable clades in target tree `t` using trees in `tref`. Only non-conflicting resolutions are kept. 
"""
function _resolve_trees_target(t::Tree, tref::Vararg{Tree}; verbose=false)

    resolved_clades = unique(vcat([sort(_resolve_trees_ref(t, tr)) for tr in tref]...))
    # Of all these clades, keep only ones that do not intersect 
    todel = Int64[]
    for (i,c1) in enumerate(resolved_clades), (j,c2) in enumerate(resolved_clades[i+1:end])
        if !isempty(intersect(c1,c2)) && intersect(c1,c2) != c1 && intersect(c1,c2) != c2
            push!(todel, i)
            push!(todel, j+i)
        end
    end
    sort!(todel); unique!(todel)
    verbose && println("$(length(todel)) clades have resolution conflicts")
    deleteat!(resolved_clades, todel)
    verbose && println("$(length(resolved_clades)) clade can be resolved")
    # 
    return resolved_clades
end

"""
    _resolve_trees_ref!(t::Tree, tref::Tree, rtau)

Return resolvable clades in  `t` using `tref`. Clades in `tref` that could be made clades in `t` by adding one internal node are resolved. 
"""
function _resolve_trees_ref!(t::Tree, tref::Tree, rtau)
    
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
            if !isclade(nodes) && is_unresolved_clade(nodes, checkclade=false) 
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

"""
    _resolve_trees_ref(t::Tree, tref::Tree)

Return resolvables clades in `t` using `tref`. Similar to `_resolve_trees_ref!` but return clades as label lists instead of changing `t`. 
"""
function _resolve_trees_ref(t::Tree, tref::Tree)
    
    !share_labels(t, tref) && @error "Can only resolve trees that share leaf nodes."
    #
    resolved_clades = Array{String,1}[]
    for (i,cl) in enumerate(keys(tref.lleaves))
        nroot = tref.lleaves[cl]
        while !nroot.isroot 
            # Go up one clade in `tref`
            nroot = nroot.anc
            nlabel = node_leavesclade_labels(nroot) # Clade in `tref`
            nodes = [t.lleaves[x] for x in nlabel]
            if !isclade(nodes) && is_unresolved_clade(nodes, checkclade=false) 
                push!(resolved_clades, sort([x.label for x in nodes]))
            else
                break
            end
        end
    end
    return unique(resolved_clades) # Note: the same clade may appear multiple times here, since we can reach it from different leaves.
end

"""
    resolve_null_branches!(t)

Branches shorter than `tau` are set to `tau`. 
"""
function resolve_null_branches!(t::Tree; tau = 1e-10, internal=false)
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
function make_unresolved_clade!(node_list::Array{<:TreeNode{T},1} ; label="", tau = 0., safe=true) where T

    if safe && !is_unresolved_clade(node_list)
        error("Input is not an unresolved clade")
    end


    mrca = lca(node_list)
    checklist = Dict((x.label)=>false for x in node_list)
    r = TreeNode(T())
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