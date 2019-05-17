export maximal_coherent_clades
# export is_coherent_clade_nodelist
# export is_coherent_clade
export name_mcc_clades!
export adjust_branchlength!
export supraMCCs
export compute_mcc_scores
export compute_mcc_scores_pairs

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


"""
Check whether a `nodelist` forms a coherent clade in all trees of `treelist`. 
i. Check that `nodelist` is a clade
ii. Find the common ancestor to `nodelist`.
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
For each clade `m` in `MCCs`: 
- Rename the root `r` of `m` to `MCC_\$(i)` or (`\$(r.label)` if `r` is a leaf) where `i` is an integer starting at `label_init`.
- Rename each non-leaf internal node of `m` to `shared_\$i_\$j` where `j` is an index specific to `m`.  

## Procedure
In an MCC internal node is defined in all trees by the clade it forms. 
"""
function name_mcc_clades!(treelist, MCCs)
    # Finding initial label
    label_init = 1
    for t in treelist
        for n in values(t.nodes)
            if match(r"MCC", n.label)!=nothing && parse(Int64, n.label[5:end]) >= label_init
                label_init = parse(Int64, n.label[5:end]) + 1
            end
        end
    end
    # println("label_init = $label_init")

    for (i,m) in enumerate(MCCs)
        cl = i + label_init - 1
        # Renaming root
        for t in treelist
            r = lca([t.lnodes[x] for x in m])
            old_label = r.label
            new_label = r.isleaf ? "$(old_label)" : "MCC_$(cl)"
            r.label = new_label
            delete!(t.lnodes, old_label)
            t.lnodes[new_label] = r
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
end

"""
"""
function adjust_branchlength!(treelist, tref, MCCs)
    # Checking that MCCs make sense before adjusting
    for (i,m) in enumerate(MCCs)
        nlist = [tref.lleaves[x] for x in m]
        if !is_coherent_clade_nodelist(nlist, (treelist..., tref))
            println(i,m)
            error("Input MCCs are not consistent clade in all trees.")
        end
    end

    # Adjusting branch lenght
    for m in MCCs
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
    supraMCCs(treelist, MCC)

Find supra MCCs: clades that are common to all trees in `treelist` and contain as few MCCs as possible (i.e. they should be direct ancestors to MCCs ideally)
Method: For each `m` in MCC 
1. Start with the root of `m` in `first(treelist)`: `r` 
2. The clade `C` defined by `m.anc` is our first candidate to a supraMCC
3. For each tree in `treelist`, check if `C` is a clade. If not, `a = mrca(C)` and `C<--clade(a)`. 
4. Iterate 3. until `C` is a clade for all trees in `treelist`. `C` is the supraMCC corresponding to `m`
"""
function supraMCCs(treelist, MCC)
    supra = Array{SupraMCC,1}(undef, 0)
    for m in MCC
        sup = SupraMCC()
        flag = true
        r = lca([first(treelist).lnodes[x] for x in m]).anc # Ancestor of `m` in one of the trees
        # if isnothing(r)
        #     r = first(treelist).root
        # end
        llist = node_leavesclade_labels(r)
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
        sup.labels = Set(llist)
        if !in(llist, [x.labels for x in supra])
            push!(supra, sup)
        end
    end
    # Composition in terms of MCCs
    for s in supra
        for m in MCC
            if prod([in(x,s.labels) for x in m])
                push!(s.MCC, lca(first(treelist).lnodes[x] for x in m).label)
            end
        end
    end
    return supra
end

"""
    compute_mcc_scores_pairs(segtrees, jointtree, MCC ; nmax = 15)

Function to score pairs of MCCs based on their removal. Outputs two scores for each pair `(m1,m2)` in MCC: 
1. average size of remaining MCCs after removing `m1` and `m2`
2. number of remaining MCCs after removing `m1` and `m2`
"""
function compute_mcc_scores_pairs(segtrees, jointtree, MCC ; nmax = 15)
    MCC_scores = Dict{Tuple{String,String}, Array{Float64,1}}()
    if length(MCC)<nmax
        for i in 1:length(MCC)
            # print("$(i)/$(length(MCC))             \r")
            for j in (i+1):length(MCC)
                m1 = MCC[i]
                m2 = MCC[j]
                jpruned = prunenodes(jointtree,cat(m1,m2,dims=1))
                tpruned = Dict(k=>prunenodes(segtrees[k],cat(m1,m2,dims=1)) for k in keys(segtrees))
                for (k,t) in tpruned
                    tpruned[k] = remove_internal_singletons(t)
                end
                Iterating._resolve_trees!(tpruned, jpruned)
                MCCn = maximal_coherent_clades(collect(values(tpruned)))
                li = lca(first(values(segtrees)).lnodes[x] for x in m1).label
                lj = lca(first(values(segtrees)).lnodes[x] for x in m2).label
                MCC_scores[li,lj] = [mean([length(x) for x in MCCn]), length(MCCn)]
            end
        end
    end
    return MCC_scores
end

"""
    compute_mcc_scores(segtrees, jointtree, MCC ; nmax = 50)

Function to score MCCs based on the effect of their removal. Outputs two scores for `m` in MCC: 
1. average size of remaining MCCs after removing `m`
2. number of remaining MCCs after removing `m`
"""
function compute_mcc_scores(segtrees, jointtree, MCC ; nmax = 100)
    MCC_scores = Dict{String, Array{Float64,1}}()
    if length(MCC)<nmax
        for (i,m) in enumerate(MCC)
            print("$(i)/$(length(MCC))             \r")
            jpruned = prunenodes(jointtree,m)
            tpruned = Dict(k=>prunenodes(segtrees[k],m) for k in keys(segtrees))
            for (k,t) in tpruned
                tpruned[k] = remove_internal_singletons(t)
            end
            Iterating._resolve_trees!(tpruned, jpruned)
            MCCn = maximal_coherent_clades(collect(values(tpruned)))
            MCC_scores[lca(first(values(segtrees)).lnodes[x] for x in m).label] = [mean([length(x) for x in MCCn]), length(MCCn)]
        end
    else
        # @warn()
    end
    return MCC_scores
end







#####################################################################################
################################### OLD STUFF #######################################
#####################################################################################

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