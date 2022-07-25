

function draw_ARG(trees:: Vector{Tree{TreeTools.MiscData}}, MCCs_dict::Dict{Set{String}, Vector{Vector{String}}};
    label_nodes = true, draw_connections = false) 
    
    first_tree = trees[1]
    tree_list = trees[2:end]
    rec_sites_dict = prepare_trees!(first_tree, tree_list, MCCs_dict)
    draw_ARG(first_tree, tree_list, rec_sites_dict; label_nodes=label_nodes, draw_connections = draw_connections)
end

function prepare_trees!(first_tree::Tree{TreeTools.MiscData}, tree_list::Vector{Tree{TreeTools.MiscData}}, 
    MCCs_dict::Dict{Set{String}, Vector{Vector{String}}})

    for tree in tree_list
        TreeKnit.assign_mccs!(TreeKnit.get_mcc_map(MCCs_dict[Set([tree.label, first_tree.label])]), tree)
    end
    assign_all_mccs!(first_tree, tree_list, MCCs_dict)
    recombination_sites = TreeKnit.get_recombination_sites(first_tree, tree_list, MCCs_dict)
    return recombination_sites
end


function draw_ARG(
    first_tree::Tree{TreeTools.MiscData}, tree_list::Union{Nothing, Vector{Tree{TreeTools.MiscData}}},
    recombination_sites::Union{Nothing, Vector{Dict{Int, Tuple{TreeNode, TreeNode}}}}; label_nodes = true,
    draw_connections = false)

    @assert isnothing(recombination_sites) && isnothing(tree_list) || !isnothing(recombination_sites) && !isnothing(tree_list)
    if !isnothing(tree_list)
        @assert haskey(first_tree.root.data.dat, "mcc") && all([haskey(tree.root.data.dat, "mcc") for tree in tree_list])
    else
        tree_list = []
        recombination_sites = Dict()
    end
    gr(size=(1200,1000), xtickfontsize=13, ytickfontsize=13, xguidefontsize=16, yguidefontsize=16, legendfontsize=10, dpi=100, grid=(:y, :gray, :solid, 1, 0.4));

    tree = first_tree
    colors = [:red, :blue, :brown, :green]
    labels = ["tree "*tree.label]
    labels = append!(labels, ["tree "*tree.label for tree in tree_list])


    function get_x_positions(tree)
        depths = Dict{TreeNode,Float16}()
        for (name, node) in tree.lnodes
            if !ismissing(node.tau)
                depths[node] = distance(node, tree.root; topological=false)
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
    # The function draw_clade closes over the axes object
    # fig = plt.figure()
    # axes = fig.add_subplot(1, 1, 1)

    # Add margins around the tree to prevent overlapping the axes
    xmax = maximum(values(x_posns))
    plot(xlims=(-0.05 * xmax, 1.25 * xmax), ylims=(0.2, (maximum(values(y_posns)) + 0.8)), ylabel="taxa", xlabel="branch length", label="", legend = :bottomright, title="Dendrogram of ARG")


    function get_x_positions_recombination(recombination_sites)
        depths_recom = Dict{Int,Float16}()
        for (mcc, n) in recombination_sites
            if !ismissing(n[1].tau)
                depths_recom[mcc] = distance(n[1], tree.root; topological=false)
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
        recombination_sites_coord = Dict{Int, Tuple{Float16, Float16, String}}()
        for mcc in keys(recombination_sites[pos])
            recombination_sites_coord[mcc] = (x_recom[pos][mcc], y_posns[recombination_sites[pos][mcc][1]] -(pos+1)*epsilon, "r"*string(mcc))
        end
        append!(recombination_sites_coord_list, [recombination_sites_coord])

        for (mcc, coord) in recombination_sites_coord

            #Add a single point using marker formatted as desired:
            plot!([coord[1]],[coord[2]], markershape=:x,
                markercolor=colors[pos], markersize=5, markerstrokewidth = 3, label="")
            if label_nodes
                annotate!((coord[1],coord[2],text(coord[3], 10, :bottom, :left, colors[pos])))
            end
        end
    end

    if draw_connections

        recom_label_list = []

        for pos in 1:length(tree_list)
            recom_labels = Dict{Int, Tuple{Float16, Float16, String}}()
            for (mcc, n) in recombination_sites[pos]
                found = false
                child_list = []
                node = n[1]
                go_up = 0
                while !found
                    for child in node.anc.child
                        if child.label !=n[1].label
                            if haskey(tree.lnodes, child.label) && !occursin("RESOLVED", child.label)
                                tree_node = tree.lnodes[child.label]
                                t_n = tree_node
                                for i in 1:go_up
                                    if t_n != tree.root
                                        t_n = t_n.anc
                                    end
                                end
                                if t_n == tree.root
                                    t_p = t_n
                                else
                                    t_p = t_n.anc
                                end
                                recom_labels[mcc] = (x_posns[t_p]/2 + x_posns[t_n]/2, y_posns[t_n] -(pos+1)*epsilon, "r'"*string(mcc))
                                found = true
                                break
                            else
                                child_list+= child
                            end
                        end
                    end
                    if !found
                        node = pop!(child_list)
                        go_up += 1
                    end
                end
            end
            for mcc in keys(recombination_sites_coord_list[pos])
                plot!([recom_labels[mcc][1]], [recom_labels[mcc][2]], markershape=:x,
                markercolor=colors[pos], markersize=5, markerstrokewidth = 3, label="")
                annotate!((recom_labels[mcc][1],recom_labels[mcc][2],text(recom_labels[mcc][3], 10, :bottom, :left, colors[pos])))
                plot!([recom_labels[mcc][1], recombination_sites_coord_list[pos][mcc][1]], [recom_labels[mcc][2], recombination_sites_coord_list[pos][mcc][2]], color=colors[pos], linestyle=:dash, label="")
                
            end
        end
    end

    function draw_clade_lines(;
        orientation="horizontal", y_here=0, x_start=0, x_here=0, y_bot=0, y_top=0, color=:black, lw=.1, l=""
    )
        if orientation=="horizontal"
            display(plot!([x_start; x_here], [y_here; y_here], lw=lw, lc=color, label=l))
        else
            display(plot!([x_here; x_here], [y_bot; y_top], lw=lw, lc=color, label=l))
        end
    end

    function draw_clade(clade, x_start, lw, x_posns, y_posns; 
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
        if isroot(clade)
            if isnothing(pos)
                label_ = labels[1]
            else
                label_ = labels[pos+1]
            end
        else
            label_ = ""
        end
        # Draw a horizontal line from start to here
        draw_clade_lines(
            orientation="horizontal",
            y_here=y_here,
            x_start=x_start,
            x_here=x_here,
            color=color,
            lw=lw, l=label_
        )
        # Add node/taxon labels
        if !recombination
            label = string(clade.label)
            if !isnothing(label) && label_nodes
                annotate!((x_here,y_here,text(label, 10, :center, :left, :black)))
            end
        end
        if !isempty(clade.child)
            # Draw a vertical line connecting all children
            y_top = y_posns[clade.child[1]]
            y_bot = y_posns[clade.child[end]]
            if recombination
                y_top -= (pos+1)*epsilon
                y_bot -= (pos+1)*epsilon
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
                        draw_clade(child, x_here, lw, x_posns, y_posns, color=color, recombination=recombination, recombination_sites =recombination_sites, pos=pos)
                    end
                end
            else
                for child in clade.child
                    draw_clade(child, x_here, lw, x_posns, y_posns, color=color, recombination=recombination)
                end
            end
        end
    end

    draw_clade(tree.root, 0, 2, x_posns, y_posns)
    for pos in 1:length(tree_list)
        draw_clade(tree.root, 0, 2, x_posns, y_posns, recombination=true, pos=pos, recombination_sites =Set(keys(recombination_sites[pos])))
        for mcc in keys(recombination_sites[pos])
            draw_clade(recombination_sites[pos][mcc][1], x_recom[pos][mcc], 2, x_posns, y_posns, recombination=true, pos=pos)
        end
    end
    savefig("arg_plot.png")

end
