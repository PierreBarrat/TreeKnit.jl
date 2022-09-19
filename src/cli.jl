"""
treeknit

# Arguments

- `nwk1`: Newick file for first tree
- `nwk2`: Newick file for second tree

# Options

- `-o, --outdir <arg>`: directory for results; Example `-o=treeknit_results`
- `-g, --gamma <arg>`: value of γ; Example `-g=2`
- `--seq-lengths <arg>`: length of the sequences. Example: `--seq-length "1500 2000"`
- `--n-mcmc-it <arg>`: number of MCMC iterations per leaf; default 25

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
	n_mcmc_it::Int = 25,
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

	# Reading sequence lengths
	sl = try
		local sl = map(i->parse(Int, i), split(seq_lengths, " "))
		@assert length(sl) == 2
		sl
	catch err
		@error "Unrecognized format for `--seq-lengths`.
Should be of the form `--seq-lengths \"1500 2000\"`"
		error(err)
	end

	# Setting up OptArgs
	oa = OptArgs(;
		γ = gamma,
		likelihood_sort = !no_likelihood,
		resolve = !no_resolve,
		nMCMC = n_mcmc_it,
		seq_lengths = sl,
		verbose=true,
	)

	#
	@info "Parameters: $oa"

	@info "Inferring MCCs...\n"
	out = @timed computeMCCs(t1, t2, oa; naive)
	MCCs = out[1]
	@info "Found $(length(MCCs)) MCCs (runtime $(out[2]))\n"

	verbose && println()

	@info "Resolving trees based on found MCCs..."
	t1, t2, rS = resolve_strict(t1, t2, MCCs)
	TreeTools.ladderize!(t1)
	sort_polytomies_strict!(t1, t2, MCCs)
	@info "Resolved $(length(rS[1])) splits in $(nwk1) and $(length(rS[1])) splits in $(nwk2)\n"

	verbose && println()

	@info "Building ARG from trees and MCCs..."
	arg, rlm, lm1, lm2 = SRG.arg_from_trees(t1, t2, MCCs)
	@info "Found $(length(arg.hybrids)) reassortments in the ARG.\n"

	verbose && println()

	# Write output
	@info "Writing results in $(outdir)"
	write_mccs(outdir * "/" * "MCCs.dat", MCCs)
	out_nwk1, out_nwk2 = make_output_tree_names(nwk1, nwk2)
	write_newick(outdir * "/" * out_nwk1, t1)
	write_newick(outdir * "/" * out_nwk2, t2)
	write(outdir * "/" * "arg.nwk", arg)
	write_rlm(outdir * "/" * "nodes.dat", rlm)

	close(io)

	println()
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

function make_output_tree_names(nwk1, nwk2)
	fn = [basename(nwk) for nwk in (nwk1, nwk2)]
	name1, name2 = if fn[1] == fn[2]
		name, ext = splitext(fn[1])
		d1 = split(dirname(nwk1), '/')[end]
		d2 = split(dirname(nwk2), '/')[end]
		name1, name2 = if d1 == d2
			name * "_1.resolved" * ext, name * "_2.resolved" * ext
		else
			name * "_$(d1).resolved" * ext, name * "_$(d2).resolved" * ext
		end
		@warn "The two input trees have the same filename. Writing output as:
		$nwk1 --> $name1
		$nwk2 --> $name2"
		name1, name2
	else
		f1, ext1 = splitext(fn[1])
		f2, ext2 = splitext(fn[2])
		f1 * ".resolved" * ext1, f2 * ".resolved" * ext2
	end

	return name1, name2
end


