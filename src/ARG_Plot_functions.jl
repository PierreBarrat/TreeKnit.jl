using PyCall
using Conda
Conda.add("biopython")

function ARGPlot(tree_list_input, MCC_list; draw_connections=false, tree_names=nothing)
## use a python wrapper to plot trees using Bio.Phylo
py"""
from Bio import Phylo
from Bio import MissingPythonDependencyError
import matplotlib.pyplot as plt
from io import StringIO

def ARGPlot(tree_list_input, MCC_list, draw_connections=False, tree_names=None):

    main_tree = Phylo.read(StringIO(tree_list_input[0]), "newick")
    tree_list = []
    for i in range(1,len(tree_list_input)):
        tree_list.append(Phylo.read(StringIO(tree_list_input[i]), 'newick'))

    print("\nARGPlot: reading trees into python")

    ## add mccs as labels to trees, similar to assign.mccs in TreeTime
    print(main_tree)
    print(tree_list)
    prepare_trees(main_tree, tree_list, MCC_list)
    print("ARGPlot: pre-preparing trees")
    ## find location of recombination events in the tree
    recombination_pairs_list = get_recombination_sites(main_tree, tree_list, MCC_list)
    print("ARGPlot: finding recombination sites")
    ## draw the ARG from one tree using a modified Phylo.draw function
    draw_ARG(main_tree, tree_list, recombination_pairs_list, draw_connections=draw_connections, tree_name_list=tree_names)
    print("ARGPlot: drawing ARG \n")
    plt.savefig('ARG_plot.png')

def print_tree(tree_string, number):
    ax = plt.axes()
    tree = Phylo.read(StringIO(tree_string), "newick")
    Phylo.draw_ascii(tree)
    #plt.title("Tree "+ str(number))
    #Phylo.draw(tree, axes=ax)


def assign_mccs_first_tree(first_tree, len_tree_list, mcc_map, one_mutation=1e-4):
    # assign MCCs to leaves
    for leaf in first_tree.get_terminals():
        leaf.child_mccs = [set([mcc_map[leaf.name][pos]]) for pos in range(len_tree_list)]
        leaf.mcc = mcc_map[leaf.name]
        leaf.branch_length = max(0.5*one_mutation, leaf.branch_length)

    # reconstruct MCCs with Fitch algorithm
    for n in first_tree.get_nonterminals(order='postorder'):
        common_mccs = [set.intersection(*[c.child_mccs[pos] for c in n]) for pos in range(len_tree_list)]
        n.branch_length = max(0.5*one_mutation, n.branch_length)
        n.child_mccs = []
        for pos in range(len_tree_list):
            if len(common_mccs[pos]):
                n.child_mccs.append(common_mccs[pos])
            else:
                n.child_mccs.append(set.union(*[c.child_mccs[pos] for c in n]))

    mcc_intersection = [set.intersection(*[c.child_mccs[pos] for c in first_tree.root]) for pos in range(len_tree_list)]
    first_tree.root.mcc = []
    for pos in range(len_tree_list):
        if len(mcc_intersection[pos]):
            first_tree.root.mcc.append(list(mcc_intersection[pos])[0])
        else:
            first_tree.root.mcc.append(None)

    for n in first_tree.get_nonterminals(order='preorder'):
        if n==first_tree.root:
            continue
        else:
            n.mcc = []
            for pos in range(len_tree_list):
                if n.up.mcc[pos] in n.child_mccs[pos]: # parent MCC part of children -> that is the MCC
                    n.mcc.append(n.up.mcc[pos])
                elif len(n.child_mccs[pos])==1:  # child is an MCC
                    n.mcc.append(list(n.child_mccs[pos])[0])
                else: # no unique child MCC and no match with parent -> not part of an MCCs
                    n.mcc.append(None)

def assign_mccs(tree, pos, mcc_map, one_mutation=1e-4):

    # assign MCCs to leaves
    for leaf in tree.get_terminals():
        leaf.child_mccs = set([mcc_map[leaf.name][pos]])
        leaf.mcc = mcc_map[leaf.name][pos]
        leaf.branch_length = max(0.5*one_mutation, leaf.branch_length)

    # reconstruct MCCs with Fitch algorithm
    for n in tree.get_nonterminals(order='postorder'):
        common_mccs = set.intersection(*[c.child_mccs for c in n])
        n.branch_length = max(0.5*one_mutation, n.branch_length)
        if len(common_mccs):
            n.child_mccs = common_mccs
        else:
            n.child_mccs = set.union(*[c.child_mccs for c in n])

    mcc_intersection = set.intersection(*[c.child_mccs for c in tree.root])
    if len(mcc_intersection):
        tree.root.mcc = list(mcc_intersection)[0]
    else:
        tree.root.mcc = None

    for n in tree.get_nonterminals(order='preorder'):
        if n==tree.root:
            continue
        else:
            if n.up.mcc in n.child_mccs: # parent MCC part of children -> that is the MCC
                n.mcc = n.up.mcc
            elif len(n.child_mccs)==1:  # child is an MCC
                n.mcc = list(n.child_mccs)[0]
            else: # no unique child MCC and no match with parent -> not part of an MCCs
                n.mcc = None


def prepare_trees(first_tree, tree_list, MCCs_list):

    #add information for later iterations
    full_tree_list = tree_list.copy()
    full_tree_list.append(first_tree)
    for tree in full_tree_list:
        tree.root.branch_length = 0

        def label_parent(n):
            parent = n
            if len(n.clades)>0:
                for child in n.clades:
                    child.up = parent
                    label_parent(child)

        n = tree.root
        label_parent(n)
        depths= n.depths(unit_branch_lengths=True)
        for n in tree.find_clades():
            n.branch_length = depths[n]
    # make a lookup for the MCCs and assign to trees
    leaf_to_MCC = {}
    for MCCs in MCCs_list:
        for mi,mcc in enumerate(MCCs):
            for leaf in mcc:
                if leaf not in leaf_to_MCC:
                    leaf_to_MCC[leaf] = [mi]
                else:
                    leaf_to_MCC[leaf].append(mi)

    for pos in range(len(tree_list)):
        assign_mccs(tree_list[pos], pos, leaf_to_MCC)
    assign_mccs_first_tree(first_tree, len(tree_list), leaf_to_MCC)

def get_recombination_sites(first_tree, tree_list, MCC_lists):

    recombination_pairs_list = []
    for pos in range(len(tree_list)):
        tree = tree_list[pos]
        checked_mccs = set()
        recombination_sites = {}
        for n in tree.get_terminals():
            if n.mcc not in checked_mccs:
                while n!=tree.root:
                    if n.mcc != n.up.mcc:
                        recombination_sites[n.mcc] = n
                        break
                    n = n.up
            checked_mccs.add(n.mcc)
            if len(checked_mccs) == len(MCC_lists[pos]):
                break
        checked_mccs_first = set()
        recombination_sites_first = {}
        for n in first_tree.get_terminals():
            if n.mcc[pos] not in checked_mccs_first:
                while n!=first_tree.root:
                    if n.mcc[pos] != n.up.mcc[pos]:
                        recombination_sites_first[n.mcc[pos]] = n
                        break
                    n = n.up
            checked_mccs.add(n.mcc[pos])
            if len(checked_mccs_first) == len(MCC_lists[pos]):
                break
        recombination_pairs = {}
        for mcc in recombination_sites_first:
            if mcc in recombination_sites:
                recombination_pairs[mcc] = (recombination_sites_first[mcc], recombination_sites[mcc])
        recombination_pairs_list += [recombination_pairs]
    return recombination_pairs_list


from Bio import MissingPythonDependencyError

def draw_ARG(
        first_tree, tree_list,
        recombination_sites,
        draw_connections = False,
        tree_name_list = None,
        label_func=str,
        do_show=True,
        show_confidence=True,
        # For power users
        axes=None,
        branch_labels=None,
        label_colors=None,
        *args,
        **kwargs):
        
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        try:
            import pylab as plt
        except ImportError:
            raise MissingPythonDependencyError(
                "Install matplotlib or pylab if you want to use draw."
            ) from None

    import matplotlib.collections as mpcollections

    tree = first_tree
    colors = ['red', 'blue', 'brown', 'green']

    # Arrays that store lines for the plot of clades
    horizontal_linecollections = []
    vertical_linecollections = []

    # Options for displaying branch labels / confidence
    def conf2str(conf):
        if int(conf) == conf:
            return str(int(conf))
        return str(conf)

    if not branch_labels:
        if show_confidence:

            def format_branch_label(clade):
                try:
                    confidences = clade.confidences
                    # phyloXML supports multiple confidences
                except AttributeError:
                    pass
                else:
                    return "/".join(conf2str(cnf.value) for cnf in confidences)
                if clade.confidence is not None:
                    return conf2str(clade.confidence)
                return None

        else:

            def format_branch_label(clade):
                return None

    elif isinstance(branch_labels, dict):

        def format_branch_label(clade):
            return branch_labels.get(clade)

    else:
        if not callable(branch_labels):
            raise TypeError(
                "branch_labels must be either a dict or a callable (function)"
            )
        format_branch_label = branch_labels

    # options for displaying label colors.
    if label_colors:
        if callable(label_colors):

            def get_label_color(label):
                return label_colors(label)

        else:
            # label_colors is presumed to be a dict
            def get_label_color(label):
                return label_colors.get(label, "black")

    else:

        def get_label_color(label):
            # if label_colors is not specified, use black
            return "black"

    # Layout

    def get_x_positions(tree):
        depths = tree.depths()
        # If there are no branch lengths, assume unit branch lengths
        if not max(depths.values()):
            depths = tree.depths(unit_branch_lengths=True)
        return depths


    def get_y_positions(tree):
        maxheight = tree.count_terminals()
        # Rows are defined by the tips
        heights = {
            tip: maxheight - i for i, tip in enumerate(reversed(tree.get_terminals()))
        }


        # Internal nodes: place at midpoint of children
        def calc_row(clade):
            for subclade in clade:
                if subclade not in heights:
                    calc_row(subclade)
            # Closure over heights
            heights[clade] = (
                heights[clade.clades[0]] + heights[clade.clades[-1]]
            ) / 2.0

        if tree.root.clades:
            calc_row(tree.root)
        return heights

    x_posns = get_x_positions(tree)
    y_posns = get_y_positions(tree)
    # The function draw_clade closes over the axes object
    if axes is None:
        fig = plt.figure()
        axes = fig.add_subplot(1, 1, 1)
    elif not isinstance(axes, plt.matplotlib.axes.Axes):
        raise ValueError(f"Invalid argument for axes: {axes}")

    def get_x_positions_recombination(recombination_sites):
        depths_recom = {}
        for mcc, n in recombination_sites.items():
            depths_recom[mcc] = n[0].branch_length/2
            parent = n[0].up
            while parent !=tree.root:
                depths_recom[mcc] += parent.branch_length
                parent = parent.up
        return depths_recom

    x_recom = [get_x_positions_recombination(recombination_sites[pos]) for pos in range(len(tree_list))]
    
    epsilon = 0.05
    recombination_sites_coord_list = []

    for pos in range(len(tree_list)):
        recombination_sites_coord = {}
        for mcc in recombination_sites[pos]:
            recombination_sites_coord[mcc] = (x_recom[pos][mcc], y_posns[recombination_sites[pos][mcc][0]] -(pos+1)*epsilon, "r" + str(mcc))
        recombination_sites_coord_list.append(recombination_sites_coord)

        for mcc, coord in recombination_sites_coord.items():
            axes.scatter(coord[0], coord[1], marker="x", color=colors[pos])
            plt.text(x=coord[0]+0.01,y =coord[1]-0.1,s=coord[2], color=colors[pos])
    
    if draw_connections:

        def get_node(tree, label):
            if "RESOLVED" in label:
                return False
            l = [n for n in tree.find_clades() if label==n.name]
            if len(l) >0:
                return l[0]
            else:
                return False

        recom_label_list = []

        for pos in range(len(tree_list)):
            recom_labels = {}
            for mcc, n in recombination_sites[pos].items():
                found = False
                child_list = []
                node = n[1]
                go_up = 0
                while not found:
                    for child in node.up.clades:
                        if child.name !=n[1].name:
                            tree_node = get_node(tree, child.name)
                            if tree_node:
                                t_n = tree_node
                                for i in range(go_up):
                                    if t_n != tree.root:
                                        t_n = t_n.up
                                if t_n == tree.root:
                                    t_p = t_n
                                else:
                                    t_p = t_n.up
                                recom_labels[mcc] = (x_posns[t_p]/2 + x_posns[t_n]/2, y_posns[t_n] -(pos+1)*epsilon, "r'" + str(mcc))
                                found = True
                                break
                            else:
                                child_list+= child
                    if not found:
                        node = child_list.pop()
                        go_up += 1
            
            for mcc in recombination_sites_coord_list[pos]:
                axes.scatter(recom_labels[mcc][0], recom_labels[mcc][1], marker="x", color=colors[pos])
                plt.text(x=recom_labels[mcc][0]+0.01,y =recom_labels[mcc][1]-0.1,s=recom_labels[mcc][2], color=colors[pos])
                plt.plot([recom_labels[mcc][0], recombination_sites_coord_list[pos][mcc][0]], [recom_labels[mcc][1], recombination_sites_coord_list[pos][mcc][1]], color=colors[pos], linewidth=plt.rcParams["lines.linewidth"], linestyle='--')

    def draw_clade_lines(
        use_linecollection=False,
        orientation="horizontal",
        y_here=0,
        x_start=0,
        x_here=0,
        y_bot=0,
        y_top=0,
        color="black",
        lw=".1",
    ):
        if not use_linecollection and orientation == "horizontal":
            axes.hlines(y_here, x_start, x_here, color=color, lw=lw)
        elif use_linecollection and orientation == "horizontal":
            horizontal_linecollections.append(
                mpcollections.LineCollection(
                    [[(x_start, y_here), (x_here, y_here)]], color=color, lw=lw
                )
            )
        elif not use_linecollection and orientation == "vertical":
            axes.vlines(x_here, y_bot, y_top, color=color)
        elif use_linecollection and orientation == "vertical":
            vertical_linecollections.append(
                mpcollections.LineCollection(
                    [[(x_here, y_bot), (x_here, y_top)]], color=color, lw=lw
                )
            )

    def draw_clade(clade, x_start, color, lw, recombination=False, recombination_sites = None, pos=None, epsilon=0.05):
        epsilon = epsilon
        x_here = x_posns[clade]
        y_here = y_posns[clade]
        if recombination:
            mcc = clade.mcc
            y_here -= (pos+1)*epsilon

        if recombination:
            color = colors[pos]
        else:
            # phyloXML-only graphics annotations
            if hasattr(clade, "color") and clade.color is not None:
                color = clade.color.to_hex()
        if hasattr(clade, "width") and clade.width is not None:
            lw = clade.width * plt.rcParams["lines.linewidth"]
        # Draw a horizontal line from start to here
        draw_clade_lines(
            use_linecollection=True,
            orientation="horizontal",
            y_here=y_here,
            x_start=x_start,
            x_here=x_here,
            color=color,
            lw=lw,
        )
        # Add node/taxon labels
        if not recombination:
            label = label_func(clade)
            if label not in (None, clade.__class__.__name__):
                axes.text(
                    x_here,
                    y_here,
                    f" {label}",
                    verticalalignment="center",
                    color=get_label_color(label),
                )
        # Add label above the branch (optional)
        conf_label = format_branch_label(clade)
        if conf_label:
            axes.text(
                0.5 * (x_start + x_here),
                y_here,
                conf_label,
                fontsize="small",
                horizontalalignment="center",
            )
        if clade.clades:
            # Draw a vertical line connecting all children
            y_top = y_posns[clade.clades[0]]
            y_bot = y_posns[clade.clades[-1]]
            if recombination:
                y_top -= (pos+1)*epsilon
                y_bot -= (pos+1)*epsilon
                l = len(clade.clades)
                if clade.clades[0].mcc[pos] != mcc[pos]:
                    if mcc[pos] is None and clade.clades[0].mcc[pos] not in recombination_sites:
                        i = 0
                    else:
                        i=1
                        while i+1<l and ((mcc[pos] is not None and clade.clades[i+1].mcc[pos]!= mcc[pos])
                        or (mcc[pos] is None and clade.clades[i+1].mcc[pos] in recombination_sites)) and y_posns[clade.clades[i+1]]<=y_here:
                            i+=1
                        y_top = min(y_posns[clade.clades[i]], y_here)
                if clade.clades[-1].mcc[pos] != mcc[pos]:
                    if mcc[pos] is None and clade.clades[-1].mcc[pos] not in recombination_sites:
                        i = 1
                    else:
                        i=2
                        while i+1<(l+1) and ((mcc[pos] is not None and clade.clades[-i-1].mcc[pos]!= mcc[pos])
                        or (mcc[pos] is None and clade.clades[-i-1].mcc[pos] in recombination_sites)) and y_posns[clade.clades[-i-1]]>=y_here:
                            i+=1
                        y_bot = max(y_posns[clade.clades[-i]], y_here)
            # Only apply widths to horizontal lines, like Archaeopteryx
            draw_clade_lines(
                use_linecollection=True,
                orientation="vertical",
                x_here=x_here,
                y_bot=y_bot,
                y_top=y_top,
                color=color,
                lw=lw,
            )
            # Draw descendents
            if recombination:
                for child in clade:
                    if child.mcc[pos] == mcc[pos] or mcc[pos] is None and child.mcc[pos] not in recombination_sites:
                        draw_clade(child, x_here, color, lw, recombination=recombination, recombination_sites =recombination_sites, pos=pos)
            else:
                for child in clade:
                    draw_clade(child, x_here, color, lw, recombination=recombination)

    draw_clade(tree.root, 0, "k", plt.rcParams["lines.linewidth"])
    for pos in range(len(tree_list)):
        draw_clade(tree.root, 0, "k", plt.rcParams["lines.linewidth"], recombination=True, pos=pos, recombination_sites =set(recombination_sites[pos].keys()))
        for mcc in recombination_sites[pos]:
            draw_clade(recombination_sites[pos][mcc][0], x_recom[pos][mcc], "k", plt.rcParams["lines.linewidth"], recombination=True, pos=pos)
    # If line collections were used to create clade lines, here they are added
    # to the pyplot plot.
    for i in horizontal_linecollections:
        axes.add_collection(i)
    for i in vertical_linecollections:
        axes.add_collection(i)

    # Aesthetics

    try:
        name = tree.name
    except AttributeError:
        pass
    else:
        if name:
            axes.set_title(name)
    axes.set_xlabel("branch length")
    axes.set_ylabel("taxa")
    # Add margins around the tree to prevent overlapping the axes
    xmax = max(x_posns.values())
    axes.set_xlim(-0.05 * xmax, 1.25 * xmax)
    # Also invert the y-axis (origin at the top)
    # Add a small vertical margin, but avoid including 0 and N+1 on the y axis
    axes.set_ylim(max(y_posns.values()) + 0.8, 0.2)

    # Parse and process key word arguments as pyplot options
    for key, value in kwargs.items():
        try:
            # Check that the pyplot option input is iterable, as required
            list(value)
        except TypeError:
            raise ValueError(
                'Keyword argument "%s=%s" is not in the format '
                "pyplot_option_name=(tuple), pyplot_option_name=(tuple, dict),"
                " or pyplot_option_name=(dict) " % (key, value)
            ) from None
        if isinstance(value, dict):
            getattr(plt, str(key))(**dict(value))
        elif not (isinstance(value[0], tuple)):
            getattr(plt, str(key))(*value)
        elif isinstance(value[0], tuple):
            getattr(plt, str(key))(*value[0], **dict(value[1]))
    
    if tree_name_list is not None:
        from matplotlib.lines import Line2D

        lines = [Line2D([0], [0], marker='o', color='black', label=tree_name_list[0])]
        for pos in range(len(tree_list)):
            lines.append(Line2D([0], [0], marker='o', color=colors[pos], label=tree_name_list[1+pos]))
        axes.legend(handles=lines)
"""

    py"ARGPlot"(tree_list_input, MCC_list, draw_connections=draw_connections, tree_names=tree_names)
end