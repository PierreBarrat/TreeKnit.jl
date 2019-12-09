export maximal_coherent_clades
# export is_coherent_clade_nodelist
# export is_coherent_clade
export name_mcc_clades!
export adjust_branchlength!
export supraMCC
export prune_mcc_scores
export prune_mcc_scores_pairs, prune_mcc_scores_triplets
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
    for (cl,v) in checklist
        # print("$i/$(length(t.lleaves)) -- $cl  --     Found $(sum(values(checklist)))/$(length(t.lleaves))            \r")
        if !v # If leave not already visited
            # We're going to go up in all trees at the same time
            croot = map(x->x.lleaves[cl], treelist_) 
            clabel = [cl]
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
                        tcroot = [lca(c) for c in nclade]
                        if tcroot == croot # Singleton in the tree, or clade with a single node --> the algorithm is getting stuck on this node
                            croot = [x.anc for x in croot]
                        else
                            croot = tcroot
                        end
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
    end
    return mc_clades
end
maximal_coherent_clades(t...) = maximal_coherent_clades(collect(t))


"""
Check whether a `nodelist` forms a coherent clade in all trees of `treelist`. 
i. Check that `nodelist` is a clade
ii. Find the common ancestor to `nodelist`.
iii. Check whether this common ancestor is a coherent clade  
All members of `nodelist` should be leaves.
"""
function is_coherent_clade_nodelist(nodelist::Array{<:TreeNode,1}, treelist)
    if !mapreduce(x->x.isleaf, *, nodelist)
        error("All nodes in `nodelist` should be leaves.")
    end
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
        for n in values(t.nodes)
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
        for n in node_clade(r)
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
    supraMCC(treelist, MCC)

Find supra MCC: clades that are common to all trees in `treelist` and contain as few MCC as possible (i.e. they should be direct ancestors to MCCs ideally)
Method: For each `m` in MCC 
1. Start with the root of `m` in `first(treelist)`: `r` 
2. The clade `C` defined by `r.anc` is our first candidate to a supraMCC
3. For each tree in `treelist`, check if `C` is a clade. If not, `a = mrca(C)` and `C<--clade(a)`. 
4. Iterate 3. until `C` is a clade for all trees in `treelist`. `C` is the supraMCC corresponding to `m`
"""
function supraMCC(treelist, MCC)
    supra = Array{Array{String,1},1}(undef, length(MCC))
    for (i,m) in enumerate(MCC)
        r = lca([first(treelist).lnodes[x] for x in m]).anc # Ancestor of `m` in one of the trees
        llist = node_leavesclade_labels(r)
        flag = true
        while flag
            flag = false
            for t in treelist
                mapr = [t.lnodes[x] for x in llist]
                if !isclade(mapr)
                    flag = true
                    r = lca(mapr)
                    llist = node_leavesclade_labels(r)
                end
            end
        end
        supra[i] = llist
    end
    return supra
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
function reduce_to_mcc(tree::Tree, MCC)
    if !assert_mcc((tree,), MCC)
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












