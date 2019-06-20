module Interfacing


using TreeTools, RecombTools
using FastaIO


global augur = "/home/pierrebc/miniconda3/envs/augur/bin/augur"
global iqtree = "/home/pierrebc/miniconda3/envs/augur/bin/iqtree"
global treetime = "/home/pierrebc/miniconda3/envs/augur/bin/treetime"

function infertree(alignment::String, outfile::String, opt::Vararg{String})
	# Setting up temporary directory
	mkpath("tmp_infertree")
	run(`cp $alignment tmp_infertree`)


	# Inferring
	base = ["-redo", "-s", "tmp_infertree/$alignment"]
	for x in opt
		push!(base, x)
	end
	run(pipeline(`$iqtree $base`), stdout="tmp_infertree/log.txt")

	# Moving results
	run(`mv tmp_infertree/$alignment.treefile $outfile`)
	run(`rm -rd tmp_infertree`)
end

"""
	scalebranches(tree::Tree, fasta)

Uses iqtree to scale branch length. 
0. Create a temporary directory. Store sequences contained in `tree` to it for later recovery. 
1. Add outgroup `RecombTools.Iterating.root_strain` to tree
2. Writes `tree` to a file
3. Find relevant sequences in `fasta`, and writes them to a file
4. Call `iqtree -s msa -t tree -blscale` 
5. Read result, remove outgroup, re-store sequences in the `Tree` object`, and return it
"""
function scalebranches(_tree::Tree, fasta)
	# 0.
	mkpath("tmp")
	TreeTools.write_fasta("tmp/stored_sequences.fasta", _tree)
	# 1.
	tree = deepcopy(_tree)
	graftnode!(tree.root, TreeNode(label="$(RecombTools.Iterating.root_strain)"))
	tree = node2tree(tree.root)
	tree.lnodes["$(RecombTools.Iterating.root_strain)"].data.tau = 0.
	write_newick("tmp/btree.nwk", tree.root)
	# 2.
	FastaWriter("tmp/align.fasta") do f
		for (n,s) in fasta
			if in(n, [x.label for x in values(tree.leaves)])
				writeentry(f, n, s)
			end
		end
	end
	# 3.
	run(pipeline(`$iqtree -s tmp/align.fasta -t tmp/btree.nwk -blscale -o $(replace(RecombTools.Iterating.root_strain, "/"=>"_"))`,stdout="tmp/log.txt"))
	# 4.
	# run(pipeline(`sed 's/_/\//g' tmp/align.fasta.treefile`, stdout="tmp/otree.nwk"))
	out = node2tree(read_newick("tmp/align.fasta.treefile"))
	for n in values(out.leaves)
		# Replacing 2 first and last _ by / -- iqtree replaces all /s by _s
		# If a country ever has a _ in its name, I'll turn into the Hulk
		n.label = replace(n.label, "_"=>"/", count=2)
		n.label = replace(n.label[end:-1:1], "_"=>"/", count=1)[end:-1:1]
	end
	isempty(_tree.leaves[1].data.sequence) ? L = 10000 : L = length(_tree.leaves[1].data.sequence)
	delete_null_branches!(out.root, threshold = 1/L/3)
	out = node2tree(out.root)
	out = prunenode(out, RecombTools.Iterating.root_strain, propagate=false)
	fasta2tree!(out, "tmp/stored_sequences.fasta")
	run(`rm -rd tmp`)
	return out
end


"""
"""
function run_treetime(; verbose=false, aln="", tree="", out="")
	verbose ? v=1 : v=0
	mkpath("$out")
	run(pipeline(`$treetime ancestral --verbose $v --aln $aln --tree $tree --out $out`, stdout="$out/log.txt"))
end


"""
	concatenate_fasta(outfile::String, fasta::Vararg{String})

Concatenate alignments `fasta` in a single one. 
1. Check that sequences appear in all input alignments
2. 
"""
function concatenate_fasta(outfile::String, fasta::Vararg{String})
	# Build dict for 1.
	seqcount = Dict()
	for f in fasta
		for (n,s) in FastaReader(f)
			nn = split(n,'|')[1]
			seqcount[nn] = get(seqcount, nn, 0) + 1
		end
	end
	#  
	catseq = Dict()
	labels = Dict()
	for f in fasta
		for (n,s) in FastaReader(f)
			nn = split(n,'|')[1]
			if seqcount[nn] == length(fasta)
				catseq[nn] = get(catseq, nn, "")*s
				labels[nn] = n
			end
		end
	end
	#
	open(outfile,"w") do f
    	for k in keys(labels)
    		write(f, ">$(labels[k])\n$(catseq[k])\n")
    	end
	end
end

end