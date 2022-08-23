"""
treeknit

# Arguments

- `nwk_files`: Newick files vector

# Options

- `-o, --outdir <arg>`: directory for results; Example `-o=treeknit_results`
- `-g, --gamma <arg>`: value of γ; Example `-g=2`
- `--seq-lengths <arg>`: length of the sequences. Example: `--seq-length "1500 2000"`
- `--n-mcmc-it <arg>`: number of MCMC iterations per leaf; default 25

# Flags

- `--naive`: Naive inference (overrides `-g`).
- `--no-likelihood`: Do not use branch length likelihood test to sort between different MCCs
- `--no-resolve`: Do not attempt to resolve trees before inferring MCCs.
- `--force-consistent`: Force output MCCs to be consistent with each other by potentially further splitting MCCs
- `--parallel`: Run sequential multiTK with parallelization (only for 4 or more trees) 
- `-v, --verbose`: verbosity
"""
@main function treeknit(
	nwk_file1::AbstractString, nwk_file2::AbstractString,
	nwk_files::AbstractString...;
	# options
	outdir::AbstractString = "treeknit_results",
	gamma::Float64 = 2.,
	seq_lengths::AbstractString = join([string(i) for i in repeat([1], 2+ length(nwk_files))], " "),
	n_mcmc_it::Int = 25,
	# flags
	naive::Bool = false,
	no_likelihood::Bool = false,
	no_resolve::Bool = false,
	verbose::Bool = false,
	force_consistent::Bool = false,
	parallel::Bool = false,
)

	nwk_files = [nwk_file1, nwk_file2, nwk_files...]
	println("Treeknit: ")
	println("Input trees:")
	println(join(["$nwk \t" for nwk in nwk_files]))
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
	@info "Input Newick files:"*join(["$nwk \t" for nwk in nwk_files])*". Reading trees..."
	fn = get_tree_names(nwk_files)
	trees = [read_tree(nwk, label=f) for (nwk, f) in zip(nwk_files, fn)]
	for t_i in trees[2:end] 
		if !TreeTools.share_labels(trees[1], t_i)
			error("Trees must share leaves")
		end
	end

	# Reading sequence lengths
	sl = try
		local sl = map(i->parse(Int, i), split(seq_lengths, " "))
		@assert length(sl) == length(trees)
		sl
	catch err
		@error "Unrecognized format for `--seq-lengths`.
Should be of the form `--seq-lengths \"1500 2000\"`"
		error(err)
	end

	##only parallelize if there are 4 or more trees
	if length(trees)<=3
		parallel = false
	end

	# Setting up OptArgs
	oa = OptArgs(;
		γ = gamma,
		likelihood_sort = !no_likelihood,
		resolve = !no_resolve,
		nMCMC = n_mcmc_it,
		seq_lengths = sl,
		verbose=true,
		force_consist=force_consistent,
		parallel=parallel,
	)

	#
	@info "Parameters: $oa"

	@info "Inferring MCCs...\n"
	out = @timed get_infered_MCC_pairs!(trees, oa)
	MCCs = out[1]
	l = [length(m) for (key,m) in MCCs.mccs]
	@info "Found $l MCCs (runtime $(out[2]))\n"
	verbose && println()

	if length(trees) ==2
		t1 = trees[1]
		t2 = trees[2]
		@info "Building ARG from trees and MCCs..."
		arg, rlm, lm1, lm2 = SRG.arg_from_trees(t1, t2, get(MCCs, t1.label, t2.label))
		@info "Found $(length(arg.hybrids)) reassortments in the ARG.\n"
		trees[1] = t1
		trees[2] = t2
	else
		##write the ARG for 2+ trees
	end

	verbose && println()
	# Write output
	@info "Writing results in $(outdir)"
	write_mccs(outdir * "/" * "MCCs.json", MCCs)
	out_nwk = make_output_tree_names(fn)
	for i in 1:MCCs.no_trees
		write_newick(outdir * "/" * out_nwk[i], trees[i])
	end
	if length(trees)==2
		write(outdir * "/" * "arg.nwk", arg)
		write_rlm(outdir * "/" * "nodes.dat", rlm)
	end

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

function make_output_tree_names(nwk_names)

	f = [splitext(n)[1] for n in nwk_names]
	ext = [splitext(n)[2] for n in nwk_names]
	names = [f_i * "_resolved" * ext_i for (f_i, ext_i) in zip(f, ext)]

	return names
end

function get_tree_names(nwk_files)

	fn = [basename(nwk) for nwk in nwk_files]
	if unique(fn) != fn
		name, ext = splitext(fn[1])
		d = [split(dirname(nwk), '/')[end] for nwk in nwk_files]
		fn = [name * "_$(d_i)"* ext for d_i in d]
	end
	@assert unique(fn) == fn "Input trees must be identifiable by file name"
	return fn
end

