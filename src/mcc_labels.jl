using Plots

"""
get_mcc_map(MCCs::Vector{Vector{String}}, get_cluster_no =true)

Returns a dictionary of which MCC each leaf is in, if `get_cluster_no = true` returns a list of 
which MCCs contain more than one leaf. 
"""
function get_mcc_map(MCCs::Vector{Vector{String}}; get_cluster_no =false)
    if get_cluster_no
        mcc_map = Dict{String, Int}()
        cluster_no = Int[]
        for (i,mcc) in enumerate(MCCs)
            if length(mcc)>1
                append!(cluster_no, i)
            end
            for node in mcc
                mcc_map[node] = i
            end
        end
        return mcc_map, cluster_no
    else
        mcc_map = Dict{String, Int}()
        for (i,mcc) in enumerate(MCCs)
            for node in mcc
                mcc_map[node] = i
            end
        end
        return mcc_map
    end
end

function get_mcc_map(MCCs::Vector{Vector{Vector{String}}})
    sets = [union([Set([m... ]) for m in mcc]...) for mcc in MCCs]
    @assert union(sets...) == sets[1] ## make sure labels are the same in all trees
    mcc_map = Dict{String, Vector{Int}}()
    for MCC in MCCs
        for (i,mcc) in enumerate(MCC)
            for node in mcc
                if haskey(mcc_map, node)
                    append!(mcc_map[node], i)
                else
                    mcc_map[node] = [i]
                end
            end
        end
    end
    return mcc_map
end


"""
PRT!(n::TreeNode)

Assign `mcc`s to branches (i.e. their child node) by a Pre-order traversal starting at the root node `n`.
"""
function PRT!(n::TreeNode, k::Int)
	
    if isroot(n)
        n.data["mcc"] = []
        for pos in 1:k
            if !isempty(n.data["child_mccs"][pos])
                append!(n.data["mcc"], pop!(n.data["child_mccs"][pos]))
            else
                append!(n.data["mcc"], nothing)
            end
        end
    else
        n.data["mcc"] =[]
        for pos in 1:k
            if n.anc.data["mcc"][pos] in n.data["child_mccs"][pos] # parent MCC part of children -> that is the MCC
                append!(n.data["mcc"], n.anc.data["mcc"][pos])
            elseif length(n.data["child_mccs"][pos])==1  # child is an MCC
                append!(n.data["mcc"], pop!(n.data["child_mccs"][pos]))
            else # no unique child MCC and no match with parent -> not part of an MCCs
                append!(n.data["mcc"], nothing)
            end
        end
    end

	delete!(n.data.dat, "child_mccs")

	if !isempty(n.child)
		for c in n.child
			PRT!(c, k)
		end
	end


end

function PRT!(t::Tree, k::Int)
	PRT!(t.root, k)
end

function PRT!(n::TreeNode)
    if isroot(n)
        if !isempty(n.data["child_mccs"])
            n.data["mcc"] = pop!(n.data["child_mccs"])
        else
            n.data["mcc"] = nothing
        end
    else
        if n.anc.data["mcc"] in n.data["child_mccs"] # parent MCC part of children -> that is the MCC
            n.data["mcc"] = n.anc.data["mcc"]
        elseif length(n.data["child_mccs"])==1  # child is an MCC
            n.data["mcc"] = pop!(n.data["child_mccs"])
        else # no unique child MCC and no match with parent -> not part of an MCCs
            n.data["mcc"] = nothing
        end
    end
    delete!(n.data.dat, "child_mccs")

    if !isempty(n.child)
        for c in n.child
            PRT!(c)
        end
    end

end

function PRT!(t::Tree)
	PRT!(t.root)
end

"""
add_mask!(filter::Union{Nothing, Vector{Vector{String}}}, t::Vararg{Tree})

Add a `mask` parameter to the tree, branches with a `mask` cannot have a recombination event occuring 
on them as they connect clades that should be together according to the input constraints (`filter`).
The filter should be in the form of an MCC, where if nodes are in the same clade this means they cannot
have a recombination event happen between them.

The function proceeds by allocating each node to the MCC it should be in using the Fitch algorithm, 
then branches which are in a MCC with 2 or more nodes are marked with `mask`.
"""
function add_mask!(filter::Union{Nothing, Vector{Vector{String}}}, t::Vararg{Tree})
	
	if isnothing(filter)
		return []
	end
	
	mcc_map, cluster_no = get_mcc_map(filter, get_cluster_no =true)
	# assign MCCs to leaves
    for tree in t
	    assign_mccs!(mcc_map, tree)

		for n in POT(tree)
			if n.data["mcc"] in cluster_no && (isroot(n) || n.data["mcc"]== n.anc.data["mcc"])
				n.data["mask"] = true
			else
				n.data["mask"] = false
			end
		end
	end
end

function assign_mccs!(mcc_map::Dict{String, Int}, t::Vector{Tree{TreeTools.MiscData}}) 
	
	# assign MCCs to leaves
	for tree in t
		for leaf in tree.lleaves
			leaf.second.data["child_mccs"] = Set([mcc_map[leaf.second.label]])
			leaf.second.data["mcc"] = mcc_map[leaf.second.label]
		end

		# reconstruct MCCs with Fitch algorithm
		for n in POT(tree)
			if !n.isleaf
				common_mccs = intersect([c.data["child_mccs"] for c in n.child]...)
				if !isempty(common_mccs)
					n.data["child_mccs"] = common_mccs
				else
					n.data["child_mccs"] = union([c.data["child_mccs"] for c in n.child]...)
				end
			end
		end

		PRT!(tree)

	end
end

function assign_all_mccs!(tree, len_tree_list, mcc_map)
    # assign MCCs to leaves
    for leaf in tree.lleaves
        leaf.second.data["child_mccs"] = [Set([mcc_map[leaf.second.label][pos]]) for pos in 1:len_tree_list]
        leaf.second.data["mcc"] = mcc_map[leaf.second.label]
    end

    # reconstruct MCCs with Fitch algorithm
    for n in POT(tree)
        if !n.isleaf
            common_mccs = [intersect([c.data["child_mccs"][pos] for c in n.child]...) for pos in 1:len_tree_list]
            n.data["child_mccs"] = Set{Int}[]
            for pos in 1:len_tree_list
                if !isempty(common_mccs[pos])
                    append!(n.data["child_mccs"], [common_mccs[pos]])
                else
                    append!(n.data["child_mccs"], [union([c.data["child_mccs"][pos] for c in n.child]...)])
                end
            end
        end
    end

	PRT!(tree, len_tree_list)
    
end

function get_recombination_sites(first_tree, tree_list, MCC_lists)
    
    len_trees = length(tree_list)
    recombination_pairs_list = []
    for pos in 1:len_trees
        tree = tree_list[pos]
        checked_mccs = Set{Int}()
        recombination_sites = Dict{Int, TreeNode}()
        for (key, n) in tree.lleaves
            if n.data["mcc"] ∉ checked_mccs
                while n!=tree.root
                    if n.data["mcc"] != n.anc.data["mcc"]
                        recombination_sites[n.data["mcc"]] = n
                        break
                    end
                    n = n.anc
                end
            end
            push!(checked_mccs, n.data["mcc"])
            if length(checked_mccs) == length(MCC_lists[pos])
                break
            end
        end
        checked_mccs_first = Set{Int}()
        recombination_sites_first = Dict{Int, TreeNode}()
        for (key, n) in first_tree.lleaves
            if n.data["mcc"][pos] ∉ checked_mccs_first
                while n!=first_tree.root
                    if n.data["mcc"][pos] != n.anc.data["mcc"][pos]
                        recombination_sites_first[n.data["mcc"][pos]] = n
                        break
                    end
                    n = n.anc
                end
            end
            push!(checked_mccs, n.data["mcc"][pos])
            if length(checked_mccs_first) == length(MCC_lists[pos])
                break
            end
        end
        recombination_pairs = Dict{Int, Tuple{TreeNode, TreeNode}}()
        for mcc in keys(recombination_sites_first)
            if mcc in keys(recombination_sites)
                recombination_pairs[mcc] = (recombination_sites_first[mcc], recombination_sites[mcc])
            end
        append!(recombination_pairs_list, [recombination_pairs])
        end
    end
    return recombination_pairs_list
end



function draw_ARG(
        first_tree, tree_list,
        recombination_sites;
        draw_connections = false,
        tree_name_list = nothing)

    tree = first_tree
    colors = [:red, :blue, :brown, :green]


    function get_x_positions(tree)
        depths = Dict{TreeNode,Float16}()
        for (name, node) in tree.lnodes
            if !ismissing(node.tau)
                distance(node, tree.root; topological=false)
            else
                depths[node] = node_depth(node)
            end
        end
        return depths
    end


    function get_y_positions(tree)
        maxheight = length(keys(tree.lleaves))
        # Rows are defined by the tips
        heights = Dict{TreeNode,Float16}()
        for (i, tip) in enumerate(POTleaves(tree))
            heights[tip] = i
        end

        # Internal nodes: place at midpoint of children
        function calc_row(clade)
            for c in clade.child
                if !haskey(heights, c)
                    calc_row(c)
                end
            end
            # Closure over heights
            heights[clade] = (
                heights[clade.child[1]] + heights[clade.child[end]]
            ) / 2.0
        end
        if !isempty(tree.root.child)
            calc_row(tree.root)
        end
        return heights
    end

    x_posns = get_x_positions(tree)
    y_posns = get_y_positions(tree)
    print(y_posns)
    # The function draw_clade closes over the axes object
    # fig = plt.figure()
    # axes = fig.add_subplot(1, 1, 1)

    # Add margins around the tree to prevent overlapping the axes
    xmax = maximum(values(x_posns))
    plot(xlims=(-0.05 * xmax, 1.25 * xmax), ylims=(0.2, (maximum(values(y_posns)) + 0.8)), ylabel="taxa", xlabel="branch length")


    function get_x_positions_recombination(recombination_sites)
        depths_recom = Dict{Int,Float16}()
        for (mcc, n) in recombination_sites
            if !ismissing(n[1].tau)
                depths_recom[mcc] = distance(n[1], n[1].root; topological=false)
            else
                depths_recom[mcc] = node_depth(n[1])
            end
            if n[1] !=tree.root
                depths_recom[mcc] -= distance(n[1], n[1].anc; topological=ismissing(n[1].tau))/2
            end
        end
        return depths_recom
    end

    x_recom = [get_x_positions_recombination(recombination_sites[pos]) for pos in 1:length(tree_list)]
    
    epsilon = 0.05
    recombination_sites_coord_list = []

    for pos in 1:length(tree_list)
        recombination_sites_coord = Dict()
        for mcc in keys(recombination_sites[pos])
            recombination_sites_coord[mcc] = (x_recom[pos][mcc], y_posns[recombination_sites[pos][mcc][1]] -(pos+1)*epsilon, "r"*string(mcc))
        end
        append!(recombination_sites_coord_list, recombination_sites_coord)

        for (mcc, coord) in recombination_sites_coord

            #Add a single point using marker formatted as desired:
            plot!([coord[1]],[coord[2]], markershape=:x,
                markercolor=colors[pos], markersize=10)
            annotate!((coord[1] + epsilon,coord[2],text(coord[3], 14, :bottom, :left, colors[pos])))
        end
    end
    
    # if draw_connections

    #     def get_node(tree, label):
    #         if "RESOLVED" in label:
    #             return False
    #         l = [n for n in tree.find_clades() if label==n.name]
    #         if len(l) >0:
    #             return l[0]
    #         else:
    #             return False

    #     recom_label_list = []

    #     for pos in range(len(tree_list)):
    #         recom_labels = {}
    #         for mcc, n in recombination_sites[pos].items():
    #             found = False
    #             child_list = []
    #             node = n[1]
    #             go_up = 0
    #             while not found:
    #                 for child in node.up.clades:
    #                     if child.name !=n[1].name:
    #                         tree_node = get_node(tree, child.name)
    #                         if tree_node:
    #                             t_n = tree_node
    #                             for i in range(go_up):
    #                                 if t_n != tree.root:
    #                                     t_n = t_n.up
    #                             if t_n == tree.root:
    #                                 t_p = t_n
    #                             else:
    #                                 t_p = t_n.up
    #                             recom_labels[mcc] = (x_posns[t_p]/2 + x_posns[t_n]/2, y_posns[t_n] -(pos+1)*epsilon, "r'" + str(mcc))
    #                             found = True
    #                             break
    #                         else:
    #                             child_list+= child
    #                 if not found:
    #                     node = child_list.pop()
    #                     go_up += 1
            
    #         for mcc in recombination_sites_coord_list[pos]:
    #             axes.scatter(recom_labels[mcc][0], recom_labels[mcc][1], marker="x", color=colors[pos])
    #             plt.text(x=recom_labels[mcc][0]+0.01,y =recom_labels[mcc][1]-0.1,s=recom_labels[mcc][2], color=colors[pos])
    #             plt.plot([recom_labels[mcc][0], recombination_sites_coord_list[pos][mcc][0]], [recom_labels[mcc][1], recombination_sites_coord_list[pos][mcc][1]], color=colors[pos], linewidth=plt.rcParams["lines.linewidth"], linestyle='--')

    function draw_clade_lines(;
        orientation="horizontal", y_here=0, x_start=0, x_here=0, y_bot=0, y_top=0, color=:black, lw=.1
    )
        if orientation=="horizontal"
            display(plot!([x_start; x_here], [y_here; y_here], lw=lw, lc=color, legend=false))
        else
            print("print y location")
            print(y_bot)
            print(y_top)
            display(plot!([x_here; x_here], [y_bot; y_top], lw=lw, lc=color, legend=false))
        end
    end

    function draw_clade(clade, x_start, lw; 
        color= :black, recombination=false, recombination_sites = nothing, pos=nothing, epsilon=0.05
    )
        epsilon = epsilon
        x_here = x_posns[clade]
        y_here = y_posns[clade]
        if recombination
            mcc = clade.data["mcc"]
            y_here -= (pos+1)*epsilon
        end

        if recombination
            color = colors[pos]
        end
        # Draw a horizontal line from start to here
        draw_clade_lines(
            orientation="horizontal",
            y_here=y_here,
            x_start=x_start,
            x_here=x_here,
            color=color,
            lw=lw,
        )
        # Add node/taxon labels
        if !recombination
            label = string(clade.label)
            if !isnothing(label)
                annotate!((x_here,y_here,text(label, 14, :center, :left, :black)))
            end
        end
        if !isempty(clade.child)
            # Draw a vertical line connecting all children
            y_top = y_posns[clade.child[1]]
            y_bot = y_posns[clade.child[end]]
            if recombination
                y_top -= (pos+1)*epsilon
                y_bot -= (pos+1)*epsilon
                print("print initial values")
                print(y_top)
                print(y_bot)
                print("\n")
                l = length(clade.child)
                if clade.child[1].data["mcc"][pos] != mcc[pos]
                    if isnothing(mcc[pos]) && clade.child[1].data["mcc"][pos] ∉ recombination_sites
                        i = 1
                    else
                        i=2
                        while i<l && ((!isnothing(mcc[pos]) && clade.child[i+1].data["mcc"][pos]!= mcc[pos])
                        || (isnothing(mcc[pos]) && clade.child[i+1].data["mcc"][pos] ∈ recombination_sites)) && y_posns[clade.child[i+1]]<= y_here
                            i+=1
                        end
                        print("i")
                        print(i)
                        y_top = min(y_posns[clade.child[i]], y_here)
                    end
                end
                if clade.child[end].data["mcc"][pos] != mcc[pos]
                    if isnothing(mcc[pos]) && clade.child[end].data["mcc"][pos] ∉ recombination_sites
                        i = 0
                    else
                        i=1
                        while (i+1)<(l-1) && ((!isnothing(mcc[pos]) && clade.child[end-i].data["mcc"][pos]!= mcc[pos])
                        || (isnothing(mcc[pos]) && clade.child[end-i].data["mcc"][pos] ∈ recombination_sites)) && y_posns[clade.child[end-i]]>=y_here
                            i+=1
                        end
                        y_bot = max(y_posns[clade.child[end-i]], y_here)
                    end
                end
            end
            # Only apply widths to horizontal lines, like Archaeopteryx
            draw_clade_lines(
                orientation="vertical",
                x_here=x_here,
                y_bot=y_bot,
                y_top=y_top,
                color=color,
                lw=lw,
            )
            # Draw descendents
            if recombination
                for child in clade.child
                    if child.data["mcc"][pos] == mcc[pos] || isnothing(mcc[pos]) && child.data["mcc"][pos] ∉ recombination_sites
                        draw_clade(child, x_here, lw, color=color, recombination=recombination, recombination_sites =recombination_sites, pos=pos)
                    end
                end
            else
                for child in clade.child
                    draw_clade(child, x_here, lw, color=color, recombination=recombination)
                end
            end
        end
    end

    draw_clade(tree.root, 0, 2)
    for pos in 1:length(tree_list)
        draw_clade(tree.root, 0, 2, recombination=true, pos=pos, recombination_sites =Set(keys(recombination_sites[pos])))
        for mcc in keys(recombination_sites[pos])
            draw_clade(recombination_sites[pos][mcc][1], x_recom[pos][mcc], 2, recombination=true, pos=pos)
        end
    end
    savefig("arg_plot.png")
end


    
    # if tree_name_list is not None:
    #     from matplotlib.lines import Line2D

    #     lines = [Line2D([0], [0], marker='o', color='black', label=tree_name_list[0])]
    #     for pos in range(len(tree_list)):
    #         lines.append(Line2D([0], [0], marker='o', color=colors[pos], label=tree_name_list[1+pos]))
    #     axes.legend(handles=lines)