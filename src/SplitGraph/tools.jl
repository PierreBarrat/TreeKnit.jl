export trees2graph

"""
	trees2graph(t::Vararg{TreeTools.Tree})

Create a `Graph` object `g` from a collection of trees. All trees should have the same leaf nodes.
1. An array of `LeafNode` from the leaves of trees
2. Add internal nodes of each tree `t` to `g` calling `tree2graph!(g,t,k)`, where `k` is the position of `t` in the list
"""
function trees2graph(t::Vararg{TreeTools.Tree})
	treelist = t
	# checklabels(treelist)


	# labels = collect(keys(first(treelist).lleaves))
	# labels_to_int = Dict(k=>findfirst(==(k), labels) for k in labels)
	K = length(treelist)
	labels = collect(keys(first(treelist).lleaves))
	labels_to_int = Dict(l=>k for (k,l) in enumerate(labels))
	leaves = [
		LeafNode(;
			index=i,
			conf=Int[i],
			anc=Array{SplitNode,1}(undef,K)
		) for i in 1:length(labels)
	]
	lleaves = Dict(label => leaf for (leaf, label) in zip(leaves, labels))
	g = Graph(;
		labels=labels,
		labels_to_int=labels_to_int,
		leaves=leaves,
		lleaves=lleaves,
		K=K,
	)
	for (k,t) in enumerate(treelist)
		tree2graph!(g, t, k) # Adding tree t to color k of graph
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
function tree2graph!(g::Graph, t::TreeTools.Tree, k::Int)
	if length(t.lleaves) != length(g.leaves)
		@error("`t.lleaves` and `g.leaves` do not have the same length.")
	end
	tree2splitnodes!(g, t.root, SplitNode(), k)
end


"""
	tree2splitnodes!(g::Graph, r::TreeTools.TreeNode, sr::SplitNode, k::Int)

Add children of `TreeNode` `r` to internal nodes of `g`. The `SplitNode` object `sr` corresponding to `r` should exist, except if `r` is root.
Recursively call `tree2splitnodes!` on children of `r`.
"""
function tree2splitnodes!(g::Graph, r::TreeTools.TreeNode, sr::SplitNode, k::Int)
	# If r is root, create a custom `sr`. Otherwise, use input one.
	if r.isroot
		sr = SplitNode(;
			anc = nothing,
			child = Array{GraphNode,1}(undef, length(r.child)),
			color = k,
			conf = treenode2conf(g,r),
			isroot = true,
		)
		push!(g.internals, sr)
	end
	# Special case of single leaf tree
	if r.isleaf
		push!(sr.child, g.lleaves[r.label])
		g.lleaves[r.label].anc[k] = sr
	end
	#
	for (i,c) in enumerate(r.child)
		if !c.isleaf
			sc = SplitNode(;
				anc = sr,
				child = Array{GraphNode,1}(undef, length(c.child)),
				color = k,
				conf = treenode2conf(g,c),
				isroot = false
			)
			tree2splitnodes!(g, c, sc, k)
			sr.child[i] = sc
			push!(g.internals, sc)
		else
			sr.child[i] = g.lleaves[c.label]
			if isassigned(g.lleaves[c.label].anc, k)
				error("Leaf $(c.label) has ancestor of color $k already assigned")
			end
			g.lleaves[c.label].anc[k] = sr
		end
	end
end



"""
COULD BE MADE FASTER WITH A MESSAGE PASSING ALG
"""
treenode2conf(g::Graph, n) = Int[g.labels_to_int[x.label] for x in POTleaves(n)]


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
