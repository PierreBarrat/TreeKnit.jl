export maximal_coherent_clades
export name_mcc_clades!
export adjust_branchlength!
export reduce_to_mcc

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
- Increase MCC by adding children of multiforcations one by one. I wish to keep this code as basic as possible: it should just find regions of perfect topologic compatibility in all trees of `treelist`. The rest can use another function, maybe a `combine_mcc` one. 
"""
function maximal_coherent_clades(treelist)

    # Checking that trees have the same label for leaf nodes
    sh = mapreduce(t->share_labels(t[1],t[2]), *, zip(treelist[1:end-1], treelist[2:end]))
    !sh && error("Can only be used on trees that share leaf nodes.")

    # List of already visited nodes
    # treelist_ = deepcopy(treelist)
    S = Tuple(SplitList(t) for t in treelist)
    tref = treelist[1]
    checklist = Dict(k=>false for k in keys(tref.lleaves))

    # Explore one leaf at a time
    mc_clades = []
    for (cl,v) in checklist
        if !v # If leave not already visited
            # We're going to go up in all trees at the same time
            croot = [t.lleaves[cl] for t in treelist] # Root of current maximal clade, in all trees 
            clabel = [cl]
            # Initial individual, always a common clade in all trees since it's a leaf. 
            flag = true
            while flag && prod(!x.isroot for x in croot)
                nroot = [x.anc for x in croot] # Ancestors of current maximal clade in all trees
                # Each element of `nroot` defines a new set of labels corresponding to one tree. There are two possibilites 
                # (i) Those sets of labels match. In this case, we have a potential consistent clade. To check further, call `is_coherent_clade`. 
                # (ii) Otherwise, the topology of trees in `treelist` is inconsistent above `croot`. `croot` is an MCC, break.  
                # nlabel = [Set(x.label for x in POTleaves(r)) for r in nroot] # List of sets of labels
                # if prod(nlabel[i]==nlabel[1] for i in 2:length(nlabel)) # case (i)
                if mapreduce(i->S[1].splitmap[nroot[1].label] == S[i].splitmap[nroot[i].label], *, 2:length(nroot))
                    # --> `r \in nroot` is the same split in all trees
                    # nclade = [node_leavesclade(r) for r in nroot] ## See if I can do with splits
                    # if prod(is_coherent_clade_nodelist(c, treelist) for c in nclade)
                    if is_coherent_clade(nroot,S) # check if children of `r` are also same splits in all trees
                        # tcroot = [lca(c) for c in nclade]
                        if nroot == croot # Singleton in the tree, or clade with a single node --> the algorithm is getting stuck on this node
                            croot = [x.anc for x in croot]
                        else
                            croot = nroot
                        end
                        # clabel = nlabel[1]
                        clabel = S[1].leaves[S[1].splitmap[nroot[1].label].dat]
                    else
                        flag = false
                    end
                else
                    flag = false
                end
            end
            # 
            # clabel = map(x->x.label, cclade)
            map(x->checklist[x]=true, [c for c in clabel])
            push!(mc_clades, sort([c for c in clabel]))
        end
    end
    return sort(mc_clades, lt = clt)
end
maximal_coherent_clades(t...) = maximal_coherent_clades(collect(t))
function clt(x,y)
    if length(x) < length(y)
        return true
    elseif length(x) > length(y)
        return false
    else
        return x[1] < y[1]
    end
end
function sort_mccs(mccs)
    return sort([sort(x) for x in mccs], lt=clt)
end



"""
    is_coherent_clade(r::TreeNode, S::Tuple{SplitList})

Do all children of `r` correspond to the same splits?

In practice, check that the splits `S[i].splitmap[c]` for `i` in `1:length(S)` and `c` in children of `r` are all the same. 
"""
function is_coherent_clade(roots::Array{<:TreeNode,1}, S::Tuple{<:SplitList,n} where n)
    #
    if length(roots) != length(S)
        @error "`roots` and `S` do not have the same length."
    end

    # All `r` in `roots` should at least have the same number of children
    if !mapreduce(ic-> length(roots[ic].child) == length(roots[1].child), *, 2:length(roots))
        return false
    end

    # Checking clade consistency
    for cref in roots[1].child
        if cref.isleaf
            for i in 2:length(roots)
                found = false
                for c in roots[i].child
                    if c.label == cref.label
                        found = true
                        break
                    end
                end
                if !found 
                    return false
                end
            end
        else
            children = [cref]
            sref = S[1].splitmap[cref.label]
            # Looking at child `cref` for `roots[1]`
            # For every `r` in roots, there has to be a child `c` with the same split
            for i in 2:length(roots)
                for c in roots[i].child
                    if !c.isleaf && S[i].splitmap[c.label] == sref
                        push!(children, c)
                        break
                    end
                end
                if length(children) != i # nothing found
                    return false
                end
            end

            # Clade below all children found 
            if !is_coherent_clade(children, S)
                return false
            end
        end
    end
    return true
end



"""
Check whether a `nodelist` forms a coherent clade in all trees of `treelist`. 
i. Check that `nodelist` is a clade
ii. Find the common ancestor to `nodelist`.
iii. Check whether this common ancestor is a coherent clade  
All members of `nodelist` should be leaves.
"""
function is_coherent_clade_nodelist(nodelist::Array{<:TreeNode,1}, treelist)
    if length(nodelist)==1
        return is_coherent_clade(nodelist[1], treelist)
    end
    
    if isclade(nodelist)
        A = lca(nodelist)
    else
        return false
    end
    return is_coherent_clade(A, treelist)
end


"""
Check whether the clade defined by `node` is a coherent clade in all trees of `treelist`:  
- it is a clade in all trees of `treelist`
- it is a clade at all levels, *i.e.* for all children `c` of `node`, `is_coherent_clade(c, treelist)` is true
"""
function is_coherent_clade(node::TreeNode, treelist)
    # If it's a leaf, it's a coherent clade
    if node.isleaf
        return true
    end
    # Is `node` a common clade to all trees? 
    cl = map(x->x.label, POTleaves(node))
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
    name_mcc_clades!(treelist, MCC)
    
For each clade `m` in `MCC`: 
- Rename the root `r` of `m` to `MCC_\$(i)` or (`\$(r.label)` if `r` is a leaf) where `i` is an integer starting at `label_init`.
- Rename each non-leaf internal node of `m` to `shared_\$i_\$j` where `j` is an index specific to `m`.  

## Procedure
In an MCC internal node is defined in all trees by the clade it forms. 
"""
function name_mcc_clades!(treelist, MCC)
    # Finding initial label
    label_init = 1
    for t in treelist
        for n in values(t.lnodes)
            if match(r"MCC", n.label)!=nothing && parse(Int64, n.label[5:end]) >= label_init
                label_init = parse(Int64, n.label[5:end]) + 1
            end
        end
    end

    nd = Dict()
    for (i,m) in enumerate(MCC)
        cl = i + label_init - 1
        # Renaming root
        for t in treelist
            r = lca([t.lnodes[x] for x in m])
            old_label = r.label
            new_label = r.isleaf ? "$(old_label)" : "MCC_$(cl)"
            r.label = new_label
            delete!(t.lnodes, old_label)
            t.lnodes[new_label] = r
            nd[new_label] = m
        end

        # Renaming internal nodes - Using the first element of treelist to iterate through internal nodes
        r1 = lca([treelist[1].lnodes[x] for x in m])
        j = 1
        for n in node_clade(r1) 
            if n!=r1 && !n.isleaf
                # Relevant internal node. Rename it in all trees
                # `llist` acts as a common identifier for `n` in all trees
                llist = [x.label for x in node_leavesclade(n)]
                for t in treelist
                    ln = lca([t.lnodes[x] for x in llist])
                    old_label = ln.label
                    new_label = "shared_$(cl)_$j"
                    ln.label = new_label
                    delete!(t.lnodes, old_label)
                    t.lnodes[new_label] = ln
                end
                j += 1
            end
        end
    end
    return nd
end

"""
"""
function adjust_branchlength!(treelist, tref, MCC)
    # Checking that MCC make sense before adjusting
    if !assert_mcc(treelist, MCC)
        error("MCC are not common to all trees\n")
    end

    # Adjusting branch length
    for m in MCC
        r = lca([tref.lnodes[x] for x in m])
        for n in POT(r)
            if n != r
                llist = [x.label for x in node_leavesclade(n)]
                for t in treelist
                    ln = lca([t.lnodes[x] for x in llist])
                    ln.data.tau = n.data.tau 
                end
            end
        end
    end
end

"""
    assert_mcc(treelist, MCC)

Asserts whether all elements of `MCC` are consistent clades for all trees of `treelist`. Print warning if not. Return `Bool`. 
"""
function assert_mcc(treelist, MCC)
    flag = true
    for (i,m) in enumerate(MCC)
        nlist = [first(treelist).lleaves[x] for x in m]
        if !is_coherent_clade_nodelist(nlist, treelist)
            @warn "MCC $i not common to all tree: \n $m \n"
            flag = false
        end
    end
    return flag
end

"""
   reduce_to_mcc(strainlist, MCC)

Decomposes `strainlist` into mccs that compose it. 
Errors if not possible 
"""
function reduce_to_mcc(strainlist, MCC)
    mcclist = []
    for m in MCC
        tmp = intersect(m,strainlist)
        if !isempty(tmp)
            if tmp != m
                error("`strainlist` not decomposable into mccs\n")
            else
                push!(mcclist, m)
            end
        end
    end

    return mcclist
end


"""
    reduce_to_mcc(tree, MCC)

Reduce `tree` to its MCC. Returns a tree with `length(MCC)` leaves. 
"""
function reduce_to_mcc(tree::Tree, MCC ; safe=false)
    if safe && !assert_mcc((tree,), MCC)
        error("MCC are not consistent with tree.")
    end
    #
    out = deepcopy(tree)
    for m in MCC
        r = lca([out.lnodes[x] for x in m])
        if r.isroot
            return node2tree(TreeNode(r.data, isleaf=true, isroot = true, label=r.label))
        elseif !r.isleaf
            rn = TreeNode(r.data, isleaf=true, isroot = true, label=r.label)
            a = r.anc
            prunenode!(r)
            graftnode!(a, rn)
        end
    end
    return node2tree(out.root)
end

"""
    is_branch_in_mcc(n::TreeNode, mccs::Dict)
    is_branch_in_mcc(n::TreeNode, mccs::Array)


Is the branch from `n` to `n.anc` in an MCC?   
The clade defined by `n` has to intersect with an MCC, and this intersection should be strictly smaller than the mcc itself.

# Note
This can be proven. The MCC found is unique, and `n` belongs to it. 
"""
function is_branch_in_mccs(n::TreeNode, mccs::Dict)
    cl = TreeTools.node_leavesclade_labels(n)
    for mcc in values(mccs)
        if is_branch_in_mcc(n, mcc)
        # if !isempty(intersect(cl, mcc)) && !isempty(setdiff(mcc, intersect(cl, mcc)))
            return true
        end
    end
    return false
end
function is_branch_in_mccs(n::TreeNode, mccs)
    cl = TreeTools.node_leavesclade_labels(n)
    for mcc in mccs
        if is_branch_in_mcc(n, mcc)
            # if !isempty(intersect(cl, mcc)) && !isempty(setdiff(mcc, intersect(cl, mcc)))
            return true
        end
    end
    return false
end

"""
    is_branch_in_mcc(n::TreeNode, mcc::Array{<:AbstractString})

Is the branch from `n` to `n.anc` in `mcc`?  
The clade defined by `n` has to intersect with `mcc`, and this intersection should be strictly smaller `mcc`.
"""
function is_branch_in_mcc(n::TreeNode, mcc::Array{<:AbstractString,1})
    cl = TreeTools.node_leavesclade_labels(n)
    # Check if intersection is empty
    flag = false
    for n in cl 
        if in(n, mcc)
            flag = true
            break
        end
    end
    !flag && return false

    # Check that the intersection is smaller than `mcc`, i.e. `mcc` should have one leaf that cl does not have
    for n in mcc
        if !in(n, cl)
            return true
        end
    end
    return false
end
"""
    find_mcc_with_branch(n::TreeNode, mccs::Dict)

Find the mcc to which the branch from `n` to `n.anc` belongs. If `mccs` is an array, return the pair `(index, value)`. If it is a dictionary, return the pair `(key, value)`. If no such mcc exists, return `nothing`.   
Based on the same idea that `is_branch_in_mcc`. 
"""
function find_mcc_with_branch(n::TreeNode, mccs::Dict)
    cl = TreeTools.node_leavesclade_labels(n)
    for (key,mcc) in mccs
        if !isempty(intersect(cl, mcc)) && !isempty(setdiff(mcc, intersect(cl, mcc)))
            return (key, mcc)
        end
    end
    return nothing
end
function find_mcc_with_branch(n::TreeNode, mccs::Array)
    cl = TreeTools.node_leavesclade_labels(n)
    for (i,mcc) in enumerate(mccs)
        if !isempty(intersect(cl, mcc)) && !isempty(setdiff(mcc, intersect(cl, mcc)))
            return (i,mcc)
        end
    end
    return nothing
end

"""
    is_linked_pair(n1, n2, mccs)
    is_linked_pair(n1::T, n2::T, mccs::Dict{Any,Array{T,1}}) where T
    is_linked_pair(n1::T, n2::T, mccs::Array{Array{T,1},1}) where T

Can I join `n1` and `n2` through common branches only? Equivalent to: is there an `m` in `mccs` such that `in(n1,m) && in(n2,m)`? 
"""
function is_linked_pair(n1::T, n2::T, mccs::Dict{Any,Array{T,1}}) where T
    for mcc in values(mccs)
        if in(n1, mcc)
            return in(n2, mcc)
        elseif in(n2, mcc)
            return false
        end
    end
    return false
end
function is_linked_pair(n1, n2, mccs)
    for mcc in values(mccs)
        if in(n1, mcc)
            return in(n2, mcc)
        elseif in(n2, mcc)
            return false
        end
    end
    return false
end
function is_linked_pair(n1::T, n2::T, mccs::Array{Array{T,1},1}) where T
    for mcc in mccs
        if in(n1, mcc)
            return in(n2, mcc)
        elseif in(n2, mcc)
            return false
        end
    end
    return false
end

"""
    find_mcc_with_node(n::String, mccs::Array{Array{T,1},1}) where T

Find MCC to which `n` belongs.
"""
function find_mcc_with_node(n::String, mccs::Array{Array{T,1},1}) where T
    for m in mccs
        if in(n, m)
            return m
        end
    end
    return nothing
end
find_mcc_with_node(n::TreeNode, mccs) = find_mcc_with_node(n1.label, mccs)
find_mcc_with_node(n::ARGNode, mccs) = find_mcc_with_node(n1.label, mccs)

"""
Core function for computing splits in MCCs. Include root of MCC if it's not root of tree.
"""
function _splits_in_mcc(m::Array{<:AbstractString}, t::Tree, leaves::Array{<:AbstractString}, leafmap::Dict)
    # Mask corresponding to leaves in the MCC
    mask = zeros(Bool, length(leaves))
    for l in m
        mask[leafmap[l]] = true
    end
    # 
    r = lca(t.lleaves[x] for x in m)
    S = SplitList(r, leaves, mask)
    TreeTools.clean!(S, clean_root=false)
    return S
end
"""
    splits_in_mcc(m::Array{<:AbstractString}, t::Tree)

Compute splits in MCC `m` in tree `t`. Return a `SplitList` object, with `mask` corresponding to leaves in `m`.
"""
function splits_in_mcc(m::Array{<:AbstractString}, t::Tree)
    leaves = sort(collect(keys(t.lleaves)))
    leafmap = Dict(leaf=>i for (i,leaf) in enumerate(leaves))
    return _splits_in_mcc(m, t, leaves, leafmap)
end
splits_in_mcc(m::Array{<:AbstractString}, trees::Vararg{Tree}) = Tuple(splits_in_mcc(m,t) for t in trees)

"""
    splits_in_mccs(MCCs, t::Vararg{Tree})

Compute splits in `MCCs` in tree `t`. Return an array of `SplitList` objects, with `mask` corresponding to leaves in each MCC.
When called with multiple trees, return a tuple of arrays of `SplitList` objects, each of which corresponds to one of the trees (order conserved). 
"""
function splits_in_mccs(MCCs, t::Tree)
    leaves = sort(collect(keys(t.lleaves)))
    leafmap = Dict(leaf=>i for (i,leaf) in enumerate(leaves))
    return [_splits_in_mcc(m, t, leaves, leafmap) for m in MCCs]
end
splits_in_mccs(MCCs, trees::Vararg{Tree}) = Tuple(splits_in_mccs(MCCs,t) for t in trees)








