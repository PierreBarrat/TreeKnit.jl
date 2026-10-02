using Test
using TreeTools
using TreeKnit
using TreeKnit.SplitGraph
using JSON3

@testset "OptArgs(K; kwargs...) keeps all keyword arguments" begin
	oa = OptArgs(2; γ = 5., nMCMC = 7, likelihood_sort = false, parallel = true, seq_lengths = [1000, 1500])
	@test oa.γ == 5.
	@test oa.nMCMC == 7
	@test !oa.likelihood_sort
	@test oa.parallel
	@test oa.seq_lengths == [1000, 1500]
	@test oa.resolve && oa.rounds == 1 # better_MCCs defaults for K = 2

	oa = OptArgs(3; Tmax = 0.5, nT = 10)
	@test oa.seq_lengths == [1, 1, 1]
	@test !oa.resolve && oa.rounds == 1 # better_trees defaults for K > 2
	@test oa.Trange[1] ≈ 0.5 && length(oa.Trange) == 11

	oa = OptArgs(3; method = :better_MCCs, rounds = 3)
	@test oa.resolve && oa.final_no_resolve && oa.rounds == 3

	oa = TreeKnit.set_up_optargs(2, 5.0, false, 7, [1000, 1500], false, 1, false, false, false, false, false, false)
	@test oa.γ == 5. && oa.nMCMC == 7 && !oa.likelihood_sort && oa.seq_lengths == [1000, 1500]
end

@testset "naive_mccs on more than two trees" begin
	t1 = parse_newick_string("(((A,B),C),(D,E));")
	t2 = parse_newick_string("(((A,B),C),(D,E));")
	t3 = parse_newick_string("(((A,B),D),(C,E));")
	@test naive_mccs(t1, t2, t3) == [["C"], ["D"], ["E"], ["A", "B"]]
	@test naive_mccs(t1, t2) == [["A", "B", "C", "D", "E"]]
end

@testset "iter_shared" begin
	M = MCC_set(3, ["a", "b", "c"], [[["A"], ["B"]], [["A", "B"]], [["B"], ["A"]]])
	@test TreeKnit.iter_shared(M, "a") == [get(M, "a", "b"), get(M, "a", "c")]
end

@testset "inferARG" begin
	t1 = parse_newick_string("((A:1,B:1):1,(C:1,(D:1,X:1):1):1);"; label = "t1")
	t2 = parse_newick_string("((A:1,(B:1,X:1):1):1,(C:1,D:1):1);"; label = "t2")
	arg = inferARG(t1, t2)
	@test length(arg.hybrids) == 1
	@test length(leaves(t1)) == 5 # inputs are not modified
end

@testset "missing branch lengths do not decide between configurations" begin
	@test SplitGraph.branch_likelihood(missing, 1., 1, 1) == 0.
	@test SplitGraph.branch_likelihood(1., missing, 1, 1) == 0.
	# Three equally good solutions; without branch lengths, each should be found.
	t1 = parse_newick_string("((A,B,C,D),E);")
	t2 = parse_newick_string("(((A,B),C),(D,E));")
	found = Set()
	for rep in 1:60
		M = run_treeknit!(copy(t1), copy(t2), OptArgs())
		push!(found, first(values(M.mccs)))
	end
	@test in([["D"], ["A", "B", "C", "E"]], found)
end

@testset "parallel and standard MultiTreeKnit" begin
	nwk = ["((A,(B,C)),(D,E));", "((A,B,C,D),E);", "((A,B),((C,D),E));"]
	for method in (:better_trees, :better_MCCs)
		trees = [parse_newick_string(s; label = "t$i") for (i, s) in enumerate(nwk)]
		oa = OptArgs(3; method, parallel = true)
		oa_copy = deepcopy(oa)
		M = run_treeknit!(trees, oa)
		@test oa.resolve == oa_copy.resolve && oa.seq_lengths == oa_copy.seq_lengths # not modified
		for (i, j) in [(1, 2), (1, 3), (2, 3)]
			m = get(M, i, j)
			@test sort(vcat(m...)) == ["A", "B", "C", "D", "E"]
		end
	end
end

@testset "CLI: options take effect, ARG trees match nodes.dat" begin
	dir = mktempdir()
	f1, f2 = joinpath(dir, "a.nwk"), joinpath(dir, "b.nwk")
	write(f1, "((A:1,B:1):1,(C:1,(D:1,X:1):1):1);")
	write(f2, "((A:1,(B:1,X:1):1):1,(C:1,D:1):1);")
	out = joinpath(dir, "out")
	TreeKnit.treeknit(f1, f2; outdir = out, gamma = 3.0, n_mcmc_it = 11, verbosity_level = -1)
	params = JSON3.read(read(joinpath(out, "parameters.json"), String))
	@test params[:γ] == 3.0 && params[:nMCMC] == 11
	labels = Set{String}()
	for f in ("a_liberal_resolved.nwk", "b_liberal_resolved.nwk")
		union!(labels, map(label, nodes(read_tree(joinpath(out, "ARG", f)))))
	end
	for line in eachline(joinpath(out, "ARG", "nodes.dat"))
		for x in split(line, ",")[2:3]
			x == " " || @test in(x, labels)
		end
	end
end
