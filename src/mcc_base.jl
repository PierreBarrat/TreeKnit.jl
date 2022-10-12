"""
    naive_mccs(treelist)

Find sets of nodes which are:
- clades in all trees of `treelist`,
- all subclades of nodes are clades in all trees of `treelist`
  (both of these properties define consistency),
- maximal: adding a node to a set results it in not being a clade in at least
  one of the trees.

All the trees of `treelist` should share the same leaf nodes.
"""
function naive_mccs(treelist)
    # Check that trees share leaves
    sh = mapreduce(t->share_labels(t[1],t[2]), *, zip(treelist[1:end-1], treelist[2:end]))
    !sh && error("Can only be used on trees that share leaf nodes.")

    # List of splits in trees
    S = Tuple(SplitList(t) for t in treelist)

    # List of already visited nodes
    tref = treelist[1]
    checklist = Dict(k=>false for k in keys(tref.lleaves))

    # Explore one leaf at a time
    mc_clades = []
    for (cl,v) in checklist
        if !v # If leave not already visited
            # Root of current maximal clade, in all trees
            croot = [t.lleaves[cl] for t in treelist]
            # Initial individual, always a common clade in all trees since it's a leaf.
            clabel = [cl]

            # We're going to go up in all trees at the same time
            flag = true
            while flag &&  mapreduce(x->!x.isroot, *, croot; init=true) #prod(!x.isroot for x in croot)
                # Ancestors of current maximal clade in all trees
                nroot = [x.anc for x in croot]
                # Each element of `nroot` defines a set of labels corresponding to a tree.
                # There are two possibilites
                # (i) Those sets of labels match.
                ## In this case, we have a potential consistent clade.
                ## To check further, call `is_coherent_clade`.
                ## This is checked using the previously computed `SplitList`
                # (ii) Else,
                ## The topology of trees is inconsistent above `croot`.
                ## `croot` is an MCC, break.
                if mapreduce(
                        i->S[1].splitmap[nroot[1].label] == S[i].splitmap[nroot[i].label],
                        *,
                        2:length(nroot)
                    )
                    # --> `r ∈ nroot` is the same split in all trees
                    # check if children of `r` are also same splits in all trees
                    if is_coherent_clade(nroot,S)
                        if nroot == croot
                            # Singleton in the tree, or clade with a single node
                            # --> the algorithm is stuck on this node
                            croot = [x.anc for x in croot]
                        else
                            croot = nroot
                        end
                        clabel = S[1].leaves[S[1].splitmap[nroot[1].label].dat]
                    else
                        flag = false
                    end
                else
                    flag = false
                end
            end
            #
            map(x->checklist[x]=true, [c for c in clabel])
            push!(mc_clades, sort([c for c in clabel]))
        end
    end
    return sort(mc_clades, lt = clt)
end
naive_mccs(t...) = naive_mccs(t)

function coherent_by_resolution(r, super_node, sub_splits_, super_splits_, mask)
    for c_r in r.child
        for c_s in super_node.child
            if !isleaf(c_r) && !isleaf(c_s) && !arecompatible(sub_splits_.splitmap[c_r.label], super_splits_.splitmap[c_s.label], mask)
                return false
            end
        end
    end
    return true
end

function same_by_resolution(S1, S2, nroot1, nroot2, t1, t2, mask1, mask2)
    s1 = S1.splitmap[nroot1.label] 
    s2 = S2.splitmap[nroot2.label]
    if TreeTools.is_sub_split(s1, s2, mask1)
        sub_split_ = s1
        super_split_ = s2
        sub_tree = t1
        sub_splits_ = S1
        super_splits_ = S2
        super_node = nroot2
    elseif TreeTools.is_sub_split(s2, s1, mask2)
        sub_split_ = s2
        super_split_ = s1
        sub_tree = t2
        sub_splits_ = S2
        super_splits_ = S1
        super_node = nroot1
    else
        return false
    end

    ##calculate lca of all nodes in super_split_ in sub_tree
    r = sub_tree.lleaves[sub_splits_.leaves[super_split_.dat[1]]]
    for l in 2:length(s2.dat)
        r = lca(r, sub_tree.lleaves[sub_splits_.leaves[super_split_.dat[l]]])
    end
    s_sub_lca = sub_splits_.splitmap[r.label]
    if TreeTools.isequal(s_sub_lca, super_split_, mask1)
        if coherent_by_resolution(r, super_node, sub_splits_, super_splits_, mask1)
            return true
        else
            return false
        end
    elseif !TreeTools.isequal(s_sub_lca, super_split_, mask1)&& arecompatible(s_sub_lca, super_split_, mask1)
        # Consider the set of splits just below r that are subsplits of super_splits_
        # If I join those, I should get exactly super_splits_
        stmp = Split(0)
        for n in r.child
            if n.isleaf
                l = findfirst(==(n.label), super_splits_.leaves)
                if in(l, super_split_.dat)
                    TreeTools.joinsplits!(stmp, Split([l]))
                end
            else
                if TreeTools.is_sub_split(sub_splits_.splitmap[n.label], super_split_, mask1)
                    TreeTools.joinsplits!(stmp, sub_splits_.splitmap[n.label])
                end
            end
        end

        if TreeTools.isequal(stmp, super_split_, mask1)
            if coherent_by_resolution(stmp, super_node, sub_splits_, super_splits_, mask1)
                return true
            else
                return false
            end
        else
            return false
        end
    else
        return false
    end
end

function naive_mccs(t1, t2, mcc; resolve=false)

    mask1 = [leaf in mcc for leaf in sort(collect(keys(t1.lleaves)))]
    S1 = TreeTools.SplitList(t1, mask1)
    mask2 = [leaf in mcc for leaf in sort(collect(keys(t2.lleaves)))]
    S2 = TreeTools.SplitList(t2, mask2)
    # List of splits in trees
    S = (S1, S2)

    mask_dict = Dict(zip(S[1].leaves, mask1))
    assign_mask!([t1, t2], mask_dict)

    # List of already visited nodes
    checklist = Dict(k=>false for k in mcc)
    mcc_roots = [lca(t1, mcc), lca(t2, mcc)]

    # Explore one leaf at a time
    mc_clades = []
    for (cl,v) in checklist
        if !v # If leave not already visited
            # Root of current maximal clade, in all trees
            croot = [t1.lleaves[cl], t2.lleaves[cl]]
            # Initial individual, always a common clade in all trees since it's a leaf.
            clabel = [cl]

            # We're going to go up in all trees at the same time
            flag = true
            while flag &&  croot[1]!= mcc_roots[1] &&  croot[2]!= mcc_roots[2]
                # Ancestors of current maximal clade in all trees
                nroot = [x.anc for x in croot]
                while croot[1]!= mcc_roots[1] && TreeKnit.SRG.is_shared_singleton(nroot[1], croot[1], S[1])
                    croot[1] = nroot[1]
                    nroot[1] = nroot[1].anc
                end
                while  croot[2]!= mcc_roots[2] && TreeKnit.SRG.is_shared_singleton(nroot[2], croot[2], S[2]) 
                    croot[2] = nroot[2]
                    nroot[2] = nroot[2].anc
                end
                # Each element of `nroot` defines a set of labels corresponding to a tree.
                # There are two possibilites
                # (i) Those sets of labels match.
                ## In this case, we have a potential consistent clade.
                ## To check further, call `is_coherent_clade`.
                ## This is checked using the previously computed `SplitList`
                # (ii) Else,
                ## The topology of trees is inconsistent above `croot`.
                ## `croot` is an MCC, break.
                topo_equal_ = false
                if mapreduce(
                        i->TreeTools.isequal(S[1].splitmap[nroot[1].label], S[i].splitmap[nroot[i].label], mask1),
                        *,
                        2:length(nroot)
                    )
                    # --> `r ∈ nroot` is the same split in all trees
                    # check if children of `r` are also same splits in all trees
                    if is_coherent_clade(nroot,S, mask1) 
                        if nroot == croot
                            # Singleton in the tree, or clade with a single node
                            # --> the algorithm is stuck on this node
                            croot = [x.anc for x in croot]
                        else
                            croot = nroot
                        end
                        clabel = S[1].leaves[S[1].splitmap[nroot[1].label].dat]
                        topo_equal_ = true
                    else
                        flag = false
                    end
                end
                if (resolve && !topo_equal_ && same_by_resolution(S1, S2, nroot[1], nroot[2], t1, t2, mask1, mask2))
                    if nroot == croot
                        # Singleton in the tree, or clade with a single node
                        # --> the algorithm is stuck on this node
                        croot = [x.anc for x in croot]
                    else
                        croot = nroot
                    end
                    if TreeTools.is_sub_split(S2.splitmap[nroot[2].label], S1.splitmap[nroot[1].label], mask2)
                        clabel = S[1].leaves[S[1].splitmap[nroot[1].label].dat]
                        flag = true
                    else
                        clabel = S[2].leaves[S[2].splitmap[nroot[2].label].dat]
                        flag = true
                    end
                else
                    if !topo_equal_
                        flag = false
                    end
                end
            end
            #
            map(x->checklist[x]=true, [c for c in clabel])
            clabel = filter(x->x in mcc, clabel)
            push!(mc_clades, sort([c for c in clabel]))
        end
    end
    return sort(mc_clades, lt = clt)
end

#= Custom order for MCCs =#
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
    is_coherent_clade(roots::Array{<:TreeNode,1}, S::Tuple{<:SplitList,n} where n)

Do all children of `r` correspond to the same splits?

In practice, check that the splits `S[i].splitmap[c]` for `i` in `1:length(S)` and `c` in children of `r` are all the same.
"""
function is_coherent_clade(roots::Array{<:TreeNode,1}, S::Tuple{<:SplitList,n} where n)
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

function assign_mask!(tree_list::Vector{Tree{TreeTools.MiscData}}, mask_dict) 
	
	for tree in tree_list
		# assign mask to leaves
        for leaf in tree.lleaves
            leaf.second.data["mask"] = mask_dict[leaf.second.label]
        end
    
        # if any child of a node is not masked this node should also not be masked
        for n in POT(tree)
            if !n.isleaf
                mask_sum = sum([c.data["mask"] for c in n.child])
                if mask_sum >0
                    n.data["mask"] = true
                else
                    n.data["mask"] = false
                end
            end
        end

	end
end

"""
    is_coherent_clade(roots::Array{<:TreeNode,1}, S::Tuple{<:SplitList,n} where n)

Do all children of `r` correspond to the same splits?

In practice, check that the splits `S[i].splitmap[c]` for `i` in `1:length(S)` and `c` in children of `r` are all the same.
"""
function is_coherent_clade(roots::Array{<:TreeNode,1}, S::Tuple{<:SplitList,n} where n, mask::Vector{Bool})
    
    # Checking clade consistency
    for cref in roots[1].child
        cref_ = cref
        if !cref_.data["mask"]
            continue
        end
        while !cref_.isleaf && sum([c.data["mask"] for c in cref_.child]) == 1
            cref_ =  cref_.child[findall( x -> x.data["mask"], cref_.child )[1]]
        end
        if cref_.isleaf && cref_.data["mask"] ##if mask==1 it is not masked
            for i in 2:length(roots)
                found = false
                for c in roots[i].child
                    c_ = c
                    if !c_.data["mask"]
                        continue
                    end
                    while !c_.isleaf && sum([ch.data["mask"] for ch in c_.child]) == 1
                        c_ =  c_.child[findall( x -> x.data["mask"], c_.child )[1]]
                    end
                    if cref_.label == c_.label
                        found = true
                        break
                    end
                end
                if !found
                    return false
                end
            end
        else
            children = [cref_]
            cref_masked = sort([x.label for x in filter(x->x.data["mask"], [n for n in POTleaves(cref_)])])
            # Looking at child `cref` for `roots[1]`
            # For every `r` in roots, there has to be a child `c` with the same split
            for i in 2:length(roots)
                for c in roots[i].child
                    c_ = c
                    if !c_.data["mask"]
                        continue
                    end
                    while !c_.isleaf && sum([ch.data["mask"] for ch in c_.child]) == 1
                        c_ =  c_.child[findall( x -> x.data["mask"], c_.child )[1]]
                    end
                    c_masked = sort([x.label for x in filter(x->x.data["mask"], [n for n in POTleaves(c_)])])
                    if cref_masked == c_masked ## a child of each contains the same non masked nodes
                        push!(children, c)
                        break
                    end
                end
                if !(length(children) == i) # nothing found
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
    name_mcc_clades!(treelist, MCC)

For each clade `m` in `MCC`:
- Rename the root `r` of `m` to `MCC_\$(i)` or (`\$(r.label)` if `r` is a leaf) where `i`
  is an integer starting at `label_init`.
- Rename each non-leaf internal node of `m` to `shared_\$i_\$j` where `j` is an index
  specific to `m`.

"""
function name_mcc_clades!(treelist, MCC)
    # Finding initial label
    label_init = 1
    for t in treelist
        for n in values(t.lnodes)
            if match(r"MCC", n.label)!=nothing && parse(Int, n.label[5:end]) >= label_init
                label_init = parse(Int, n.label[5:end]) + 1
            end
        end
    end

    nd = Dict()
    for (i,m) in enumerate(MCC)
        cl = i + label_init - 1
        # Renaming root
        for t in treelist
            r = lca(t, m)
            old_label = r.label
            new_label = r.isleaf ? "$(old_label)" : "MCC_$(cl)"
            r.label = new_label
            delete!(t.lnodes, old_label)
            t.lnodes[new_label] = r
            nd[new_label] = m
        end
    end

    return nd
end



"""
    reduce_to_mcc(tree, MCC)

Reduce `tree` to its MCC by grouping leaves. Returns a tree with `length(MCC)` leaves.
"""
function reduce_to_mcc(tree::Tree, MCC)
    out = copy(tree)
    reduce_to_mcc!(out, MCC)
    return out
end
"""
    reduce_to_mcc!(tree, MCC)

Reduce `tree` to `MCCs` by grouping leaves.
"""
function reduce_to_mcc!(tree::Tree, MCC)
    for m in MCC
        r = lca(tree, m)
        if r.isroot
            node2tree!(tree, TreeNode(
                    r.data;
                    isleaf=true, isroot = true, label=r.label, r.tau
                ))
        elseif !r.isleaf
        	for c in reverse(r.child)
        		# prunenode!(c)
        		prunesubtree!(tree, c; remove_singletons=false)
        	end
        	r.isleaf = true
        	r.isroot = false
        	tree.lleaves[r.label] = r
        end
    end
end

