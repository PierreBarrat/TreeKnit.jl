"""
treeknit

# Arguments

- `nwk1`: Newick file for first tree
- `nwk2`: Newick file for second tree

# Options

- `-o, --outdir <arg>`: directory for results; Example `-o=treeknit_results`
- `-g, --gamma <arg>`: value of γ; Example `-g=2`
- `--seq-lengths <arg>`: length of the sequences. Example: `--seq-length "1500 2000"`
- `--n-sa-it <arg>`: number of SA iterations per temperature and per leaf.

# Flags

- `--naive`: Naive inference (overrides `-g`).
- `--no-likelihood`: Do not use branch length likelihood test to sort between different MCCs
- `--no-resolve`: Do not attempt to resolve trees before inferring MCCs.
- `-v, --verbose`: verbosity
"""
@main function treeknit(
	nwk1::AbstractString, nwk2::AbstractString;
	# options
	outdir::AbstractString = "treeknit_results",
	gamma::Float64 = 2.,
	seq_lengths::AbstractString = "1 1",
	n_sa_it::Float64 = 1.,
	# flags
	naive::Bool = false,
	no_likelihood::Bool = false,
	no_resolve::Bool = false,
	verbose::Bool = false,
)

	println("Treeknit: ")
	println("Input trees: $nwk1 \t $nwk2")
	println("Results directory: $outdir")
	println("γ: $gamma")

	# Setting up directories
	mkpath(outdir)

	# Setting loggers
	loggers = []
	## File log.txt
	io = open(outdir*"/log.txt", "w+")
	file_logger = FormatLogger(io) do io, args
		println(io, args.file, ":", args.line, " [", args.level, "] ", args.message)
	end;
	push!(loggers, file_logger)
	## If verbose, stderr
	if verbose
		push!(loggers, ConsoleLogger())
	else
		push!(loggers, MinLevelLogger(ConsoleLogger(), Logging.Warn))
	end

	global_logger(TeeLogger(loggers...))


	# Reading trees
	@info "Input Newick files: $nwk1 \t $nwk2. Reading trees..."
	t1 = read_tree(nwk1)
	t2 = read_tree(nwk2)
	if !TreeTools.share_labels(t1, t2)
		error("Trees must share leaves")
	end

	# Setting up OptArgs
	sl = map(i->parse(Int, i), split(seq_lengths, " "))
	oa = OptArgs(;
		γ = gamma,
		likelihood_sort = !no_likelihood,
		resolve = !no_resolve,
		Md = 1 / n_sa_it,
		verbose=true,
	)

	#
	@info "Parameters: $oa"

	@info "Inferring MCCs...\n"
	MCCs = computeMCCs(t1, t2, oa; naive, seqlengths = sl)
	@info "Found $(length(MCCs)) MCCs\n"

	@info "Resolving trees based on found MCCs..."
	rS = resolve!(t1, t2, MCCs)
	@info "Resolved $(length(rS[1])) splits in $(nwk1) and $(length(rS[1])) splits in $(nwk2)\n"

	@info "Building ARG from trees and MCCs..."
	arg, rlm, lm1, lm2 = SRG.arg_from_trees(t1, t2, MCCs)
	@info "Found $(length(arg.hybrids)) reassortments in the ARG.\n"

	# Write output
	@info "Writing results in $(outdir)"
	write_mccs(outdir * "/" * "MCCs.dat", MCCs)
	for (nwk, t) in zip((nwk1,nwk2),(t1,t2))
		fn = basename(nwk)
		name, ext = splitext(fn)
		write_newick(outdir * "/" * name * ".resolved" * ext, t)
	end
	write(outdir * "/" * "arg.nwk", arg)
	write_rlm(outdir * "/" * "nodes.dat", rlm)

	close(io)
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


