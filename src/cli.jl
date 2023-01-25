"""
	treeknit

We suggest two methods for running TreeKnit:
- `--better-trees`: This method is recommended for most users using >2 trees. It will attempt to resolve all trees compatibly before inferring MCCs, it will then run 1 round of 
  pair-wise treeknit inference on all tree pairs individually - not further resolving trees. This technique will produce the most accurate output trees and the most homogeneous MCCs -
  this means that if nodes are inferred to be in a shared MCC they are infact shared. This method is also the fastest.
- `--better-MCCs`: This method is recommended for users using 2 trees, or users who are more interested in inferring accurate MCCs. It will also attempt to resolve all trees compatibly
  before inferring MCCs, it will then run 1 round of sequential treeknit inference on all tree pairs, further resolving trees using the MCCs inferred in the previous treeknit calls. 
  For K>2 trees it will then run one final round of pairwise TreeKnit on all tree pairs individually - not further resolving trees. This technique will produce the most accurate MCCs, 
  but the output trees will potentially have a higher rate of inaccurate splits. 

# Arguments

- `nwk_files`: Newick files (requires 2 or more trees - note an ARG currently cannot be created for more than 2 trees)

# Options

- `-o, --outdir <arg>`: directory for results; Example `-o=treeknit_results`
- `-g, --gamma <arg>`: value of γ; Example `-g=2`
- `--seq-lengths <arg>`: length of the sequences. Example: `--seq-length "1500 2000"`
- `--n-mcmc-it <arg>`: number of MCMC iterations per leaf; default 25
- `--rounds <arg>`: Number of times to run inference on input trees. If `rounds > 1` MCCs will be re-inferred using resolved trees from the last iteration. (default: `1`)
- `--verbosity-level <arg>`: set value of verbosity. Default 0. `-v` flag sets it to 1. Set to 2 for maximum verbosity (only useful for small trees). Set to -1 for no output at all

# Flags

- `--better-trees`: Use the `--better-trees` method (default for >2 trees)
- `--better-MCCs`: Use the `--better-MCCs` method (default for 2 trees)
- `--naive`: Naive inference (overrides `-g`).
- `--no-resolve`: Do not attempt to resolve trees before inferring pairwise MCCs.
- `--liberal-resolve`: Resolve output trees as much as possible using inferred MCCs, adding splits when order of coalescence and reassortment is unclear (coalescence is set at a time prior to reassortment)
- `--resolve-all-rounds`: Resolve trees before inferring pairwise MCCs in all rounds, overrides `--no-resolve` (default for 2 trees, for more than 2 trees default is to not resolve in the final round)
- `--no-pre-resolve`: Do not compatibly resolve all trees with each other before inferring MCCs (default is to pre-resolve)
- `--no-likelihood`: Do not use branch length likelihood test to sort between different MCCs
- `--parallel`: Run sequential TreeKnit with parallelization (only used for 3 or more trees)
- `-v, --verbose`: set verbosity to 1
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
	verbosity_level::Int = 0,
	# flags
	better_trees::Bool = false,
	better_MCCs::Bool = false,
	naive::Bool = false,
	no_pre_resolve::Bool = false,
	no_resolve::Bool = false,
	liberal_resolve::Bool = false,
	resolve_all_rounds::Bool = false,
	no_likelihood::Bool = false,
	verbose::Bool = false,
	parallel::Bool = false,
	auspice_view::Bool = false
)

	# Setting up directories
	mkpath(outdir)

	# Setting loggers
	if verbose && verbosity_level == 0
		verbosity_level = 1
	end
	io_logfile = open(outdir*"/log.txt", "w+")
	loggers = get_loggers(verbosity_level, io_logfile)
	old_logger = global_logger(TeeLogger(loggers...))

	# Start
	nwk_files = [nwk_file1, nwk_file2, nwk_files...]
	@info "Treeknit: "
	@info "Input trees: " * join(["$nwk  " for nwk in nwk_files])
	@info "Results directory: $outdir"
	@info "γ: $gamma"

	# Reading trees
	@logmsg LogLevel(-1) "Reading trees..."
	fn, ext = get_tree_names(nwk_files)
	trees = [read_tree(nwk, label=f) for (nwk, f) in zip(nwk_files, fn)]
	for t_i in trees[2:end] 
		if !TreeTools.share_labels(trees[1], t_i)
			error("Trees must share leaves")
		end
	end

	# Reading sequence lengths
	seq_lengths = try
		local sl = map(i->parse(Int, i), split(seq_lengths, " "))
		@assert length(sl) == length(trees)
		sl
	catch err
		@error "Unrecognized format for `--seq-lengths`.
				Should be of the form `--seq-lengths \"1500 2000\"`"
		error(err)
	end

	# Parallelize only if there are 3 or more trees
	if length(trees)<3 && parallel
		@warn "Cannot run in parallel for 2 trees"
		parallel = false
	end
	if parallel
		@info "Running in parallel"
	end

	naive==true && @info "Doing naive inference"

	# Setting up parameters
	oa = set_up_optargs(
		length(trees),
		gamma,
		!no_likelihood,
		n_mcmc_it,
		seq_lengths,
		parallel,
		rounds,
		better_trees,
		better_MCCs,
		no_pre_resolve,
		no_resolve,
		liberal_resolve,
		resolve_all_rounds
	)

	json_string = JSON3.write(oa)
	open(outdir * "/" *"parameters.json", "w") do f
		JSON3.pretty(f, json_string)
	end

	@logmsg LogLevel(-1) "Parameters: $oa\n\n"
	# Infer MCCs

	@info "Inferring MCCs...\n"
	infered_trees = [copy(t) for t in trees]
	out = @timed run_treeknit!(infered_trees, oa; naive)
	MCCs = out[1]

	l = [length(m) for (key,m) in MCCs.mccs]
	@info "Found $l MCCs (runtime $(out[2]))\n"
	
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

	if length(trees) ==2
		mkpath(outdir*"/ARG")
		rS = resolve!(trees[1], trees[2], get(MCCs, trees[1].label, trees[2].label); strict=false)
		out_nwk = make_output_tree_names(fn, ext; ARG=true)
		for i in 1:MCCs.no_trees
			write_newick(outdir * "/ARG/" * out_nwk[i], trees[i])
		end
		@logmsg LogLevel(-1) "Building ARG from trees and MCCs..."
		arg, rlm, lm1, lm2 = SRG.arg_from_trees(trees[1], trees[2], get(MCCs, trees[1].label, trees[2].label))
		@logmsg LogLevel(-1) "Found $(length(arg.hybrids)) reassortments in the ARG.\n"
		write(outdir * "/ARG/" * "arg.nwk", arg)
		write_rlm(outdir * "/ARG/" * "nodes.dat", rlm)
	end

	close(io_logfile)

	# give back logger
	global_logger(old_logger)

	return nothing
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

function set_up_optargs(
	K::Int,
	γ::Real,
	likelihood_sort::Bool,
	nMCMC::Int,
	seq_lengths::AbstractVector{Int},
	parallel::Bool,
	rounds::Int,
	better_trees::Bool,
	better_MCCs::Bool,
	no_pre_resolve::Bool,
	no_resolve::Bool,
	liberal_resolve::Bool,
	resolve_all_rounds::Bool
)
	@logmsg LogLevel(0) ""
	@logmsg LogLevel(0) "Setting up parameters of the TreeKnit run"
	if (better_trees == true) && (better_MCCs == true)
		error("Cannot use both `--better-trees` and `--better-MCCs`")
	end

	method = if better_trees
		@logmsg LogLevel(0) "Using the `--better-trees` method"
		:better_trees
	elseif better_MCCs
		@logmsg LogLevel(0) "Using the `--better-MCCs` method"
		:better_MCCs
	else
		if K > 2
			@logmsg LogLevel(0) "Using the `--better-trees` method by default for >2 trees"
			:better_trees
		else
			@logmsg LogLevel(0) "Using the `--better-MCCs` method by default for 2 trees"
			:better_MCCs
		end
	end
	oa = OptArgs(K; method, γ, likelihood_sort, nMCMC, seq_lengths, parallel)

	## check if any additional arguments have been passed
	(no_pre_resolve == true) && (oa.pre_resolve = false)
	(no_resolve == true) && (oa.resolve = false)
	(liberal_resolve == true) && (oa.strict = false)
	(resolve_all_rounds == true) && (oa.final_no_resolve = false)
	(rounds != 1) && (oa.rounds=rounds)

	if K>2 && oa.resolve==true && oa.final_no_resolve==false
		@warn "For K>2 we do not recommend resolving in the final round of MCC inference.
		You see this warning because you used the `--resolve-all-rounds` flag."
	end

	if (oa.final_no_resolve==true) && K==2
		@warn "Default set to not resolve in final round. This can be disabled using `--resolve-all-rounds`"
	end

	@logmsg LogLevel(0) ".. Will perform $(oa.rounds) rounds of TreeKnit"
	oa.pre_resolve==true && @logmsg LogLevel(0) ".. Will preresolve all trees before MCC inference"
	if oa.resolve==true
		@logmsg LogLevel(0) ".. Will resolve trees during MCC inference (change with `--no-resolve`)"
		oa.strict==false && @logmsg LogLevel(0) ".. Will resolve ambiguous splits with the most parsimonious option (flag --liberal-resolve)"
		oa.strict==true && @logmsg LogLevel(0) ".. Will Only resolve unambiguous splits (change with `--liberal-resolve`)"
		oa.final_no_resolve==true && @logmsg LogLevel(0) ".. Will not resolve trees in the final round of MCC inference (change with `--resolve-all-rounds`)"
	else
		@logmsg LogLevel(0) ".. Will not resolve trees during MCC inference (`--no-resolve` flag)"
	end

	return oa
end



function make_output_tree_names(nwk_names, ext; ARG=false)

	f = [splitext(n)[1] for n in nwk_names]
	if ARG
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

## Logging helpers

timestamp_logger(logger) = TransformerLogger(logger) do log
  merge(log, (; message = "[$(Dates.format(now(), date_format))] - $(log.message)"))
end
# verbosity(log) = hasproperty(log.group, :verbosity) ? log.group.verbosity : 0

function get_loggers(verbosity_level, io)
	loggers = []

	# File logger
	## show only relevant part of the path
	pth(f) = begin
		regex = r"TreeKnit.*$"
		m = match(regex, f)
		if isnothing(m)
			f
		else
			m.match
		end
	end

	file_logger = begin
		FormatLogger(io) do io, args
			println(io, pth(args.file), ":", args.line, " [", args.level, "] ", args.message)
		end |>
		x -> MinLevelLogger(x, LogLevel(-1)) |>
		timestamp_logger
	end
	push!(loggers, file_logger)

	## Console log
	#=
	TransformerLogger takes log and feeds it to ConsoleLogger
	ConsoleLogger only shows logs with level >= -verbosity_level
	if log.level < -verbosity_level, it is unchanged --> it will not be shown
	if log.level >= -verbosity_level, it is bumped up to Info --> will be shown as level Info
	=#
	console_logger =  begin
		ConsoleLogger(LogLevel(-verbosity_level)) |>
		x -> TransformerLogger(x) do log
		   if log.level >= LogLevel(-verbosity_level)
		       return merge(log, (;level = Logging.Info))
		   else
		       return log
		   end
		end |>
		timestamp_logger
    end
	push!(loggers, console_logger)

	return loggers
end
