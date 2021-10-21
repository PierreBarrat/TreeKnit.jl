using ARGTools
using RecombTools
using TreeTools
using Test

n = 15
ρ = 0.1
N = 10_000

f() = f(n, ρ, N)
function f(n, ρ, N)
	arg_ref = ARGTools.SimulateARG.simulate(N, ρ/N*n, n; simtype=:kingman)
	t1, t2 = ARGTools.trees_from_ARG(arg_ref)
	MCCs = ARGTools.MCCs_from_arg(arg_ref)

	arg = RecombTools.SRG.arg_from_trees(t1, t2, MCCs)[1]
	tt1, tt2 = RecombTools.SRG.trees_from_arg(arg)
	TreeTools.remove_internal_singletons!(tt1)
	TreeTools.remove_internal_singletons!(tt2)

	S1_ref = SplitList(t1)
	S1 = SplitList(tt1)
	S2_ref = SplitList(t2)
	S2 = SplitList(tt2)

	return S1 == S1_ref && S2 == S2_ref
end

@testset "Random args" begin
	@test mapreduce(i->f(), *, 1:50, init=true)
end
