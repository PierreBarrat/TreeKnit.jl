export trees2graph

"""
	coarsegrain_trees2graph(t::Vararg{TreeTools.Tree})
"""
function coarsegrain_trees2graph(t::Vararg{TreeTools.Tree})
	treelist = collect(t)
	mcc = maximal_coherent_clades(treelist)
	mcc_names = name_mcc_clades!(treelist, mcc)
	for (i,t) in enumerate(treelist)
		treelist[i] = reduce_to_mcc(t, mcc)
	end
	return trees2graph(treelist)
end

"""
	trees2graph(t::Vararg{TreeTools.Tree})

Create a `Graph` object `g` from a collection of trees. All trees should have the same leaf nodes. 
1. An array of `LeafNode` from the leaves of trees
2. Add internal nodes of each tree `t` to `g` calling `tree2graph!(g,t,k)`, where `k` is the position of `t` in the list
"""
function trees2graph(t::Vararg{TreeTools.Tree})
	treelist = collect(t)
	checklabels(treelist)


	K = length(treelist)
	labels = collect(keys(first(treelist).lleaves))
	labels_to_int = Dict(k=>findfirst(x->x==k, labels) for k in labels)
	N = length(labels)
	leaves = [LeafNode(index=i, conf=idx2conf(i,N), anc=Array{SplitNode,1}(undef,K)) for i in 1:length(labels)]
	lleaves = Dict(labels[i]=>leaves[i] for i in 1:length(labels))
	g = Graph(labels=labels, labels_to_int=labels_to_int, leaves=leaves, lleaves=lleaves, K=K)
	for (k,t) in enumerate(treelist)
		tree2graph!(g, t, k)
	end

	return g
end


"""
"""
function trees2graph(treelist)
	return trees2graph(treelist...)
end

"""
	tree2graph!(g::Graph, t::TreeTools.Tree, k::Int64)

Add internal nodes of `t` to `Graph` `g`. `g` should already have fields `leaves`, `labels` and `labels_to_int` initialized.  
Call `tree2splitnodes(g, t.root, SplitNode(), k)` to start the process. 
"""
function tree2graph!(g::Graph, t::TreeTools.Tree, k::Int64)
	length(t.lleaves) == length(g.leaves) || @error("`t.lleaves` and `g.leaves` do not have the same length. Graph was not initialized?")
	tree2splitnodes!(g, t.root, SplitNode(), k)
end


"""
	tree2splitnodes!(g::Graph, r::TreeTools.TreeNode, sr::SplitNode, k::Int64)

Add children of `TreeNode` `r` to internal nodes of `g`. The `SplitNode` object `sr` corresponding to `r` should exist, except if `r` is root.  
Recursively call `tree2splitnodes!` on children of `r`. 
"""
function tree2splitnodes!(g::Graph, r::TreeTools.TreeNode, sr::SplitNode, k::Int64)
	if r.isroot
		sr = SplitNode(anc=nothing, child=Array{GraphNode,1}(undef, length(r.child)), color=k, conf = treenode2conf(g,r), isroot=true)
	end 

	for (i,c) in enumerate(r.child)
		if !c.isleaf
			sc = SplitNode(anc=sr, child=Array{GraphNode,1}(undef, length(c.child)), color=k, conf = treenode2conf(g,c), isroot=false)
			tree2splitnodes!(g, c, sc, k)
			sr.child[i] = sc
			push!(g.internals, sc)
		else
			sr.child[i] = g.lleaves[c.label]
			isassigned(g.lleaves[c.label].anc, k) && error("Leaf $(c.label) has ancestor of color $k already assigned")
			g.lleaves[c.label].anc[k] = sr
		end
	end
end



"""
"""
function treenode2conf(g::Graph, n::TreeTools.TreeNode)
	N = length(g.labels)
	tmp = [g.labels_to_int[x] for x in TreeTools.node_leavesclade_labels(n)]
	conf = zeros(Bool,N)
	for idx in tmp
		conf[idx] = true
	end
	return conf
end


"""
"""
function checklabels(treelist)
	t1 = first(treelist)
	for t in treelist
		if !TreeTools.share_labels(t1,t)
			@error("Trees do not share leaves nodes labels")
		end
	end
	return nothing
end

"""
"""
function idx2conf(i::Int64, N::Int64)
	conf = zeros(Int64,N)
	conf[i] = 1
	return conf
end