"""
	struct OptArgs

Storing parameters for `SplitGraph.runopt` function.

### General
- `γ::Real = 2`
- `itmax::Int64 = 15`: maximal number of iterations of naive MCCs / SA cycles.
- `likelihood_sort::Bool = true`: sort equivalent configurations using likelihood test
  based on branch length.
- `resolve::Bool = true`: try to resolve trees while finding MCCs.
- `seq_lengths`: lengths of sequences that trees were built from.
  Used in likelihood calculations.
  This is initialized from other input arguments, and defaults to sequences of length one.
### Simulated annealing
- `Md::Real = 10`:  number of SA iterations (per temperature) for a tree of `n` leaves is
  `ceil(Int, n/Md)`.
- `cooling_schedule = :geometric`: type of cooling schedule `(:geometric, :linear)`
- `Tmin::Float64 = 1e-3`: minimal temperature of SA.
- `Tmax::Float64 = 1`: maximal temperature of SA.
- `αT::Float64 = 0.95`: ratio between terms in the geometric cooling.
- `dT::Float64 = 1e-2`: temperature step in the linear cooling.
### Verbosity
- `verbose::Bool=false`: first level of verbosity
- `vv::Bool = false`: second level of verbosity
"""
@with_kw struct OptArgs
	γ::Real = 2
	itmax::Int64 = 15
	likelihood_sort::Bool = true
	resolve::Bool = true
	seq_lengths::Vector{Int} = [1, 1]
	# For the annealing
	Md::Real = 1
	Tmin::Float64 = 1e-3
	Tmax::Float64 = 1.; @assert Tmax > Tmin
	αT::Float64 = 0.95
	dT::Float64 = 1e-2
	cooling_schedule = :geometric
	Trange = get_cooling_schedule(; Tmin, Tmax, dT, αT, type=cooling_schedule)
	sa_rep::Int64 = 1
	# Verbosity
	verbose::Bool = false
	vv::Bool = false
end

function get_cooling_schedule(;
	Tmax=1., Tmin = 1e-3, dT = 1e-2, αT = 0.95, type=:geometric,
)
	if type == :geometric
		return get_geometric_cooling_schedule(Tmin, Tmax, αT)
	elseif type == :linear
		return get_linear_cooling_schedule(Tmin, Tmax, dT)
	else
		error("Unknown `cooling_schedule` field: $(type). See `?OptArgs` for allowed values.")
	end
end

function get_geometric_cooling_schedule(Tmin, Tmax, αT)
	n = ceil(Int, (log(Tmin) - log(Tmax)) / log(αT))
	return [αT^i * Tmax for i in 0:n]
end

get_linear_cooling_schedule(Tmin, Tmax, dT) = return reverse(Tmin:dT:Tmax)


