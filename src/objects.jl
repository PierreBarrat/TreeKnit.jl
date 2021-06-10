"""
	struct OptArgs

Storing parameters for `SplitGraph.runopt` function.

### General
- `γ::Real = 3`
- `itmax::Int64 = 15`: Maximal number of iterations of MCC / SA cycles
- `likelihood_sort::Bool = true`: sort equivalent configurations using likelihood test
  based on branch length.
- `resolve::Bool = true`: try to resolve trees while finding MCCs.
- `seq_lengths = ones(Int64, 2)`: lengths of sequences that trees were built from.
  Used in likelihood calculations.
### Simulated annealing
- `Md::Real = 10`:  Number of SA iterations (per temperature) for a tree of `n` leaves is
  `ceil(Int, n/Md)`
- `Tmin::Float64 = 1e-3`: Minimal temperature of SA
- `Tmax::Float64 = 1`: Maximal temperature of SA
- `dT::Float64 = 1e-2`: Temperature step
### Verbosity
- `verbose::Bool=false`: first level of verbosity
- `vv::Bool = false`: second level of verbosity
### Output
- `output = :mccs`: possible values `[:mccs, :mccs_df, :all]`. If calling `computeMCCs`,
  this should always be set to `:mccs`. Other values are only relevant when calling
  `runopt`.
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
	dT::Float64 = 5e-2
	Trange = reverse(Tmin:dT:Tmax)
	sa_rep::Int64 = 1
	# Verbosity
	verbose::Bool = false
	vv::Bool = false
	# Output
	output = :mccs
end

