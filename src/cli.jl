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
- `--rounds`: Number of times to run inference on input trees. If `rounds > 1` MCCs will be re-inferred using resolved trees from the last iteration. (default: `1`)

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

	##only parallelize if there are 3 or more trees
	if length(trees)<3 && parallel
		@warn "Cannot run in parallel for less than 3 trees (got $(length(trees)))"
		parallel = false
	end
	if parallel
		@info "Running in parallel"
	end

	naive==true && println("Performing naive inference")

	# Setting up parameters

	oa = OptArgs(;
	γ = gamma,
	likelihood_sort = !no_likelihood,
	nMCMC = n_mcmc_it,
	seq_lengths = sl,
	verbose = true,
	parallel = parallel,
	)
	oa = set_up_OptArgs!(oa, length(trees), rounds, better_trees, better_MCCs, no_pre_resolve, no_resolve, liberal_resolve, resolve_all_rounds)

	json_string = JSON3.write(oa)
	
	open(outdir * "/" *"parameters.json", "w") do f
		JSON3.pretty(f, json_string)
	end

	@info "Parameters: $oa"
	# Infer MCCs

	@info "Inferring MCCs...\n"
	infered_trees = [copy(t) for t in trees]
	out = @timed run_treeknit!(infered_trees, oa; naive)
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
		out_nwk = make_output_tree_names(fn, ext; ARG=true)
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

function set_up_OptArgs!(oa::OptArgs, K::Int, rounds::Int, better_trees::Bool, better_MCCs::Bool, no_pre_resolve::Bool, no_resolve::Bool, liberal_resolve::Bool, resolve_all_rounds::Bool)

	if (better_trees == true) && (better_MCCs == true)
		error("Cannot use both `--better-trees` and `--better-MCCs`")
	end

	if (better_trees == true) || (better_MCCs == false && K>2)
		@info "Using `--better-trees` method"
		oa.pre_resolve = true
		oa.strict = true
		oa.final_no_resolve = true
		oa.rounds = 1
		oa.resolve = false	
	elseif (better_MCCs == true) || (better_trees == false && K==2)
		@info "Using `--better-MCCs` method"
		oa.pre_resolve = true
		oa.strict = true
		if K>2
			oa.rounds = 2
			oa.final_no_resolve = true
		else
			oa.rounds = 1
			oa.final_no_resolve = false
		end
		oa.resolve = true
	end

	(no_pre_resolve==true) && (oa.pre_resolve=false) #default false or input value
	(liberal_resolve==true) && (oa.strict=false) #default false or input value
	(resolve_all_rounds==true) && (oa.final_no_resolve=false) #default false or input value
	(rounds!=1) && (oa.rounds=rounds) #default 1 or input value
	(no_resolve==true) && (oa.resolve=false) #default false or input value

	##clean up
	if (oa.final_no_resolve==true) && (oa.rounds==1) && (oa.resolve==true)
		oa.resolve=false
		oa.final_no_resolve=false
	end

	if K>2 && oa.resolve==true && oa.final_no_resolve==false
		println("WARNING: For $K) we do not recommend resolving during MCC inference in the final round.")
	end

	if (oa.final_no_resolve==true)
		println("WARNING: Default set to not resolve in final round. This can be disabled using `--resolve-all-rounds`")
	end

	println("Performing $(oa.rounds) rounds of TreeKnit")
	oa.pre_resolve==true && println("Preresolving all trees before MCC inference")
	if oa.resolve==true
		println("Resolving trees during MCC inference")
		oa.strict==false && println("Resolve ambiguous splits with the most parsimonious option")
		oa.strict==true && println("Only resolving unambiguous splits")
		oa.final_no_resolve==true && println("Not resolving trees in the final round of MCC inference")
	else
		println("Not resolving trees during MCC inference")
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
