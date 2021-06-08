"""
	struct OptArgs

Storing parameters for `SplitGraph.runopt` function.

### General
- `γ::Real = 3`
- `itmax::Int64 = 15`: Maximal number of iterations of MCC / SA cycles
- `likelihood_sort::Bool = true`: sort equivalent configurations using likelihood test (based on branch length for now).
- `resolve::Bool = true`: try to resolve trees while finding MCCs.
- `seq_lengths = ones(Int64, 2)`: lengths of sequences that trees were built from
- `crossmap_resolve::Bool = false`: Use cross-mapped mutations to resolve polytomies in the tree.  Require each leaf's data to have `:selfseq` and `:cmseq` entries.
- `crossmap_prune::Bool = false`: Use cross-mappped mutations to prune MCCs preventively. The code will look at the number of suspicious mutations at `n.data.dat[:suspicious_muts][s]` where `n = trees[s]` (`s` is assumed to be an influenza segment). Require each leaf's data to have `:selfseq` and `:cmseq` entries.
- `suspmut_threshold::Int = 1`: Minimal number of suspicious mutation to prune a branch.
### Simulated annealing
- `Md::Real = 10`:  Number of SA iterations (per temperature) for a tree of `n` leaves is `ceil(Int64, n/Md)`
- `Tmin::Float64 = 1e-3`: Minimal temperature of SA
- `Tmax::Float64 = 1`: Maximal temperature of SA
- `dT::Float64 = 1e-2`: Temperature step
### Verbosity
- `verbose::Bool=false`: first level of verbosity
- `vv::Bool = false`: second level of verbosity
### Output
- `output = :mccs`: possible values `[:mccs, :mccs_df, :all]`
"""
@with_kw struct OptArgs
	γ::Real  = 3
	itmax::Int64 = 15
	likelihood_sort::Bool = true
	resolve::Bool = true
	seq_lengths = ones(Int64, 2)
	# For the annealing
	Md::Real = 10
	Tmin::Float64 = 1e-3
	Tmax::Float64 = 1.; @assert Tmax > Tmin
	dT::Float64 = 1e-2
	Trange = reverse(Tmin:dT:Tmax)
	sa_rep::Int64 = 1
	# Verbosity
	verbose::Bool = false
	vv::Bool = false
	# Output
	output = :mccs
end
getM(n,Md) = ceil(Int64, n/Md)
