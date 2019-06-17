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
0. Create a temporary directory 
1. Writes `tree` to a file
2. Find relevant sequences in `fasta`, and writes them to a file
3. Call `iqtree -s msa -t tree -blscale` 
4. Read result and returns a `Tree` object
"""
function scalebranches(tree::Tree, fasta)
	# 0.
	mkpath("tmp")
	# 1.
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
	run(pipeline(`$iqtree -s tmp/align.fasta -te tmp/btree.nwk -blfix`,stdout="tmp/log.txt"))
	# 4.
	run(pipeline(`sed 's/_/\//g' tmp/align.fasta.treefile`, stdout="tmp/otree.nwk"))
	out = node2tree(read_newick("tmp/otree.nwk"))
	run(`rm -rd tmp`)
	return out
end


"""
"""
function run_treetime(; verbose=false, aln="", tree="", out="")
	verbose ? v=1 : v=0
	run(pipeline(`$treetime ancestral --verbose $v --aln $aln --tree $tree --out $out`, stdout="$out/log.txt"))
end

end