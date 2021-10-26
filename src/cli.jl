"""
treeknit

# Arguments

- `nwk1`: Newick file for first tree
- `nwk2`: Newick file for second tree

# Options

- `-o, --outdir <arg>`: directory for results
- `-g, --gamma <arg>`: value of γ
- `--seq-lengths <arg>`: length of the sequences used to build trees
  Used for likelihood calculations.
  *e.g.* `--seq-length "1500 2000"`.
- `--n-sa-it <arg>`: number of SA iterations per temperature and per leaf.

# Flags

- `--naive`: Naive inference (overrides `-g`).
- `--no-likelihood`: Do not use likelihood test based on branch length to sort between
  different MCCs.
- `--no-resolve`: Do not attempt to resolve trees before inferring MCCs.
"""
@main function treeknit(
	nwk1::AbstractString, nwk2::AbstractString;
	# options
	outdir::AbstractString = "treeknit_results",
	gamma::Float64 = 2.,
	seq_lengths::AbstractString = "1 1",
	n_sa_it::Float64 = 0.1,
	# flags
	naive::Bool = false,
	no_likelihood::Bool = false,
	no_resolve::Bool = false,
)
	println("Treeknit: ")
	println("Input trees: $nwk1 \t $nwk2")
	println("Results directory: $outdir")
	println("γ: $gamma")

	# Setting up inputs
	t1 = read_tree(nwk1)
	t2 = read_tree(nwk2)

	sl = map(i->parse(Int, i), split(seq_lengths, " "))
	oa = OptArgs(;
		γ = gamma,
		likelihood_sort = !no_likelihood,
		resolve = !no_resolve,
		Md = 1 / n_sa_it,
	)

	# Setting up directories
	mkpath(outdir)

	#
	MCCs = computeMCCs(t1, t2, oa; naive, seqlengths = sl)
	resolve!(t1, t2, MCCs)
	arg, rlm, lm1, lm2 = SRG.arg_from_trees(t1, t2, MCCs)

	# Write output
	write_mccs(outdir * "/" * "MCCs.dat", MCCs)
	write_newick(outdir * "/" * nwk1 * "_resolved", t1)
	write_newick(outdir * "/" * nwk2 * "_resolved", t2)
	write(outdir * "/" * "arg.nwk", arg)
	write_rlm(outdir * "/" * "nodes.dat", rlm)
end

function write_rlm(filename, rlm)
	open(filename, "w") do io
		for (n, (k,v)) in enumerate(rlm)
			write(io, k)
			write(io, ",")
			for (i,vv) in enumerate(v)
				if isnothing(vv)
					write(io, " ")
				else
					write(io, vv)
				end
				if i != length(v)
					write(io, ",")
				end
			end
			if n != length(rlm)
				write(io, "\n")
			end
		end
	end
end

# @cast function tmp(;x = "1 1")
# 	xn = map(i->parse(Int, i), split(x, " "))
# 	println(xn[1])
# 	println(xn[2])
# end

# """
# ArgParse example implemented in Comonicon.

# # Arguments

# - `x`: an argument

# # Options

# - `--opt1 <arg>`: an option
# - `-o, --opt2 <arg>`: another option

# # Flags

# - `-f, --flag`: a flag
# """
# @cast function example1(x; opt1=1, opt2::Int=2, flag=false)
#     println("Parsed args:")
#     println("flag=>", flag)
#     println("arg=>", x)
#     println("opt1=>", opt1)
#     println("opt2=>", opt2)
# end

# @main
