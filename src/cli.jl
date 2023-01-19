"""
	treeknit

# Arguments

- `nwk_files`: Newick files (requires 2 or more trees - note an ARG currently cannot be created for more than 2 trees)

# Options

- `-o, --outdir <arg>`: directory for results; Example `-o=treeknit_results`
- `-g, --gamma <arg>`: value of γ; Example `-g=2`
- `--seq-lengths <arg>`: length of the sequences. Example: `--seq-length "1500 2000"`
- `--n-mcmc-it <arg>`: number of MCMC iterations per leaf; default 25
- `--rounds`: Number of times to run inference on input trees. If `rounds > 1` MCCs will be re-inferred using resolved trees from the last iteration. (default: `1`)

# Flags

- `--naive`: Naive inference (overrides `-g`).
- `--no-likelihood`: Do not use branch length likelihood test to sort between different MCCs
- `--no-resolve`: Do not attempt to resolve trees before inferring pairwise MCCs.
- `--liberal-resolve`: Resolve output trees as much as possible using inferred MCCs, adding splits when order of coalescence and reassortment is unclear (coalescence is set at a time prior to reassortment)
- `--resolve-all-rounds`: Resolve trees before inferring pairwise MCCs in all rounds, overrides `--no-resolve` (default for 2 trees, for more than 2 trees default is to not resolve in the final round)
- `--no-pre-resolve`: MultiTreeKnit flag - Do not compatibly resolve all trees with each other before inferring MCCs (default is to pre-resolve)
- `--parallel`: MultiTreeKnit flag - Run sequential MultiTreeKnit with parallelization (only used for 4 or more trees)
- `-v, --verbose`: verbosity
- `--auspice-view`: return ouput files for auspice
"""
@main function treeknit(
	nwk_file1::AbstractString, nwk_file2::AbstractString, nwk_files::AbstractString...;
	# options
	outdir::AbstractString = "treeknit_results",
	gamma::Real = OptArgs().γ,
	seq_lengths::AbstractString = join([string(i) for i in repeat([1], 2+ length(nwk_files))], " "),
	n_mcmc_it::Int = OptArgs().nMCMC,
	rounds::Int = 1,
	# flags
	naive::Bool = false,
	no_pre_resolve::Bool = false,
	no_likelihood::Bool = false,
	no_resolve::Bool = false,
	liberal_resolve::Bool = false,
	resolve_all_rounds::Bool = false,
	verbose::Bool = false,
	parallel::Bool = false,
	auspice_view::Bool = false
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
	fn, ext = get_tree_names(nwk_files)
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
	if length(trees)<=3 && parallel
		@warn "Not running in parallel for less than 4 trees (got $(length(trees)))"
		parallel = false
	end

	pre_resolve = !no_pre_resolve

	println("Performing $rounds rounds of TreeKnit")
	if resolve_all_rounds || (length(trees)==2 && !no_resolve)
		resolve_all_rounds = true
		no_resolve = false
		println("Resolving tree topology prior to inference")
	end
	if no_resolve
		println("Not resolving tree topology prior to inference")
		final_no_resolve = true 
	end
	if (length(trees)>2 && !resolve_all_rounds && !no_resolve)
		final_no_resolve = true
		println("Resolving tree topology prior to inference in all rounds but final round")
	else
		final_no_resolve = false
	end

	# Setting up OptArgs
	oa = OptArgs(;
		γ = gamma,
		likelihood_sort = !no_likelihood,
		resolve = !no_resolve,
		strict=!liberal_resolve,
		final_no_resolve = final_no_resolve,
		rounds = rounds, 
		nMCMC = n_mcmc_it,
		seq_lengths = sl,
		verbose = true,
		parallel = parallel,
		pre_resolve = pre_resolve,
	)

	json_string = JSON3.write(oa)
	
	open(outdir * "/" *"parameters.json", "w") do f
		JSON3.pretty(f, json_string)
	end

	#
	@info "Parameters: $oa"

	@info "Inferring MCCs...\n"
	infered_trees = [copy(t) for t in trees]
	out = @timed MTK.get_infered_MCC_pairs!(infered_trees, oa; naive)
	MCCs = out[1]

	l = [length(m) for (key,m) in MCCs.mccs]
	@info "Found $l MCCs (runtime $(out[2]))\n"
	verbose && println()
	
	# Write output
	@info "Writing results in $(outdir)"
	write_mccs(outdir * "/" * "MCCs.json", MCCs)
	out_nwk = make_output_tree_names(fn, ext)
	for i in 1:MCCs.no_trees
		write_newick(outdir * "/" * out_nwk[i], infered_trees[i])
	end

	if auspice_view
		write_auspice_json(outdir * "/", infered_trees, MCCs)
	end
	verbose && println()

	if length(trees) ==2
		mkpath(outdir*"/ARG")
		rS = resolve!(trees[1], trees[2], get(MCCs, trees[1].label, trees[2].label))
		out_nwk = make_output_tree_names(fn, ext; liberal=true)
		for i in 1:MCCs.no_trees
			write_newick(outdir * "/ARG/" * out_nwk[i], trees[i])
		end
		@info "Building ARG from trees and MCCs..."
		arg, rlm, lm1, lm2 = SRG.arg_from_trees(trees[1], trees[2], get(MCCs, trees[1].label, trees[2].label))
		@info "Found $(length(arg.hybrids)) reassortments in the ARG.\n"
		write(outdir * "/ARG/" * "arg.nwk", arg)
		write_rlm(outdir * "/ARG/" * "nodes.dat", rlm)
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

function make_output_tree_names(nwk_names, ext; liberal=false)

	f = [splitext(n)[1] for n in nwk_names]
	if liberal
		names = [f_i * "_liberal_resolved" * ext for f_i in f]
	else
		names = [f_i * "_resolved" * ext for f_i in f]
	end

	return names
end


function get_tree_names(nwk_files)

	fn = [basename(nwk) for nwk in nwk_files]
	name, ext = splitext(fn[1])
	fn = [splitext(f)[1] for f in fn]
	if !allunique(fn)
		d = [split(dirname(nwk), '/')[end] for nwk in nwk_files]
		fn = [name * "_$(d_i)" for d_i in d]
	end
	@assert allunique(fn) "Input trees must be identifiable by file name"
	return (fn, ext)
end
