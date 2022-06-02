function copysplitlist(r::TreeNode, leaves::Array{T,1}, leafmap::Dict) where T
    S = SplitList(
	    leaves,
		Array{Split,1}(undef,0),
		zeros(Bool, length(leaves)),
		Dict{eltype(leaves), Split}(),
	)
	copysplitlist!(S, r, leaves, leafmap)
	return S
end

function copysplitlist!(S::SplitList, r::TreeNode, leaves::Array{T,1}, leafmap::Dict) where T
	if r.label in leaves
		s = Split([leafmap[r.label]])
	else
		L = 0
		for c in r.child
			sc = copysplitlist!(S, c, leaves, leafmap)
			L += length(sc)
		end
		s = Split(L)
		i = 1
		for c in r.child
			if c.label in leaves
				s.dat[i] = leafmap[c.label]
				i += 1
			else
				sc = S.splitmap[c.label]
				TreeTools._joinsplits!(s,sc,i)
				i += length(sc)
			end
		end
		sort!(s.dat)
		unique!(s.dat)
		push!(S.splits, s)
		S.splitmap[r.label] = s
	end
	return s
end


function is_coherent_clade_copy_pair(roots::Array{<:TreeNode,1}, S::Tuple{<:SplitList,n} where n, leaves)
    # All `r` in `roots` should at least have the same number of non-masked children
    #if length([c if c.data["copy"][loc_pair[1]] for roots[2].child]) != 
    #        length([c if c.data["copy"][loc_pair[2]] for roots[1].child])
    #    return false
    #end

    # Checking clade consistency
    for cref in roots[1].child
        if cref.label in leaves
            found = false
            for c in roots[2].child
                if c.label == cref.label
                    found = true
                    break
                end
            end
            if !found
                return false
            end
        else
            children = [cref]
            sref = S[1].splitmap[cref.label]
            # Looking at child `cref` for `roots[1]`
            # For every `r` in roots, there has to be a child `c` with the same split
            for c in roots[2].child
                if c.label ∉ leaves && S[2].splitmap[c.label] == sref
                    push!(children, c)
                    break
                end
            end
            if length(children) != 2 # nothing found
                return false
            end

            # Clade below all children found
            if !is_coherent_clade_copy_pair(children, S, leaves)
                return false
            end
        end
    end
    return true
end

"""
    naive_mccs(treelist)

Find sets of nodes which are:
- clades in all trees of `treelist` / current tree copies in `treelist`
- all subclades of nodes are clades in all trees of `treelist`
  (both of these properties define consistency),
- maximal: adding a node to a set results it in not being a clade in at least
  one of the trees.

All the trees / tree copies of `treelist` should share the same leaf nodes. 
These are either specified in the dictionary `copyleaves` or defined as the shared 
leaves in  `treelist`
"""
function naive_mccs(treelist::Vector{Tree{T}}, copyleaves::Union{Nothing, Set{String}}) where T

    if isnothing(copyleaves)
        # Check that trees share leaves
        sh = mapreduce(t->share_labels(t[1],t[2]), *, zip(treelist[1:end-1], treelist[2:end]))
        !sh && error("Can only be used on trees that share leaf nodes.")
        copyleaves = treelist[1].lleaves
    end
    # List of splits in trees/ copy pairs
    leaves = sort(collect(copyleaves))
    leafmap = Dict(leaf=>i for (i,leaf) in enumerate(leaves))
    S = Tuple(copysplitlist(t.root, leaves, leafmap) for t in treelist)

    # List of already visited nodes
    checklist = Dict(k=>false for k in leaves)

    # Explore one leaf at a time
    mc_clades = []
    for (cl,v) in checklist
        if !v # If leave not already visited
            # Root of current maximal clade, in all trees
            croot = [t.lnodes[cl] for t in treelist]
            # Initial individual, always a common clade in all trees
            clabel = [cl]
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
                    if is_coherent_clade_copy_pair(nroot,S, leaves)
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

function naive_mccs(treelist::Vector{Tree{T}}) where T
    return naive_mccs(treelist, nothing)
end

function naive_mccs(t1::Tree{T}, t2::Tree{T}, tn::Vararg{Tree{T}}) where T
    return naive_mccs([t1, t2, tn...], nothing)
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

"""
    name_mcc_clades!(treelist, MCC)

For each clade `m` in `MCC`:
- Rename the root `r` of `m` to `MCC_\$(i)` or (`\$(r.label)` if `r` is a leaf) where `i`
  is an integer starting at `label_init`.
- Rename each non-leaf internal node of `m` to `shared_\$i_\$j` where `j` is an index
  specific to `m`.

"""
function name_mcc_clades!(treelist, copyleaves, MCC)
    # make label the same for MCC tree pair roots
    # Finding initial label
    label_init = 1
    for t in treelist
        for n in values(t.lnodes)
            if match(r"MCC", n.label)!=nothing && parse(Int, n.label[5:end]) >= label_init
                label_init = parse(Int, n.label[5:end]) + 1
            end
        end
    end
    
    # remove nodes inside MCCs from copy (they have been "cut off")
    j = 0
    for tree_pair in keys(copyleaves)
        next_j = 0
        nodes_to_be_removed_from_copy = []
        nodes_to_be_added_to_copy = []
        for (i,m) in enumerate(MCC[tree_pair])
            cl = i + label_init - 1 +j
            r_i = lca(treelist[tree_pair[1]], m)
            r_j = lca(treelist[tree_pair[2]], m)
            if r_i.label ∉ copyleaves[tree_pair]
                nodes_to_be_removed_from_copy_1 = remove_node_from_copy!(r_i, tree_pair[2], copyleaves[tree_pair])
                nodes_to_be_removed_from_copy_2 = remove_node_from_copy!(r_j, tree_pair[1], copyleaves[tree_pair])
                @assert nodes_to_be_removed_from_copy_1 == nodes_to_be_removed_from_copy_2
                append!(nodes_to_be_removed_from_copy, nodes_to_be_removed_from_copy_1)
                
                new_label = "MCC_$(cl)"
                next_j +=1
                r_i_old_label, r_j_old_label = r_i.label, r_j.label
                r_i.label, r_j.label = new_label, new_label
                delete!(treelist[tree_pair[1]].lnodes, r_i_old_label)
                delete!(treelist[tree_pair[2]].lnodes, r_j_old_label)
                treelist[tree_pair[1]].lnodes[new_label] = r_i
                treelist[tree_pair[2]].lnodes[new_label] = r_j
                push!(nodes_to_be_added_to_copy, new_label)
            end
        end
        j += next_j 
        # remove leaves inside an MCC from copyleaves and replace with root of MCC
        for n in nodes_to_be_removed_from_copy
            delete!(copyleaves[tree_pair], n)
        end
        for n in nodes_to_be_added_to_copy
            push!(copyleaves[tree_pair], n)
        end
    end

end

function remove_node_from_copy!(treenode, pos, copyleaves; internal_copyleaves=Set())
    if treenode.label ∉ copyleaves
        for c in treenode.child
            c.data.dat["copy"][pos] = 0
            if c.label ∉ copyleaves
                push!(internal_copyleaves, remove_node_from_copy!(treenode, pos, copyleaves; internal_copyleaves))
            else
                push!(internal_copyleaves, c.label)
            end
        end
    end
    return internal_copyleaves
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

