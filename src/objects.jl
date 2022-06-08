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
- `nMCMC::Int = 25`: The number *total* number of MCMC steps (swaps) for a tree of `n` leaves
	is `nMCMC*n`. The number of MCMC steps at one temperature is `nMCMC * n / nT`.
- `cooling_schedule = :geometric`: type of cooling schedule `(:geometric, :linear, :acos)`
- `Tmin::Float64 = 0.05`: minimal temperature of SA.
- `Tmax::Float64 = 0.8`: maximal temperature of SA.
- `nT::Int = 3000`: number of steps in the cooling schedule
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
	nMCMC::Int = 50
	Tmin::Float64 = 0.05; @assert Tmin > 0
	Tmax::Float64 = 1; @assert Tmax > Tmin
	nT::Int = 100
	cooling_schedule = :geometric
	Trange = get_cooling_schedule(Tmin, Tmax, nT, type=cooling_schedule)
	sa_rep::Int64 = 1
	# Verbosity
	verbose::Bool = false
	vv::Bool = false
end

function get_cooling_schedule(Tmin, Tmax, nT; type=:geometric)
	if type == :geometric
		return get_geometric_cooling_schedule(Tmin, Tmax, nT)
	elseif type == :linear
		return get_linear_cooling_schedule(Tmin, Tmax, nT)
	elseif type == :acos
		return get_acos_cooling_schedule(Tmin, Tmax, nT)
	else
		error("Unknown `cooling_schedule` field: $(type). See `?OptArgs` for allowed values.")
	end
end

function get_geometric_cooling_schedule(Tmin, Tmax, nT)
	α = exp((log(Tmin) - log(Tmax)) / nT)
	n = ceil(Int, (log(Tmin) - log(Tmax)) / log(α))
	return [α^i * Tmax for i in 0:n]
end

function get_acos_cooling_schedule(Tmin, Tmax, nT)
	f(x, K=1.5) = if x < 0.5
		0.5 + 2^(K-1)*abs(acos(2*x-1)/3.14-0.5)^K
	elseif x == 0.5
		0.5
	else
		0.5 - 2^(K-1)*abs(acos(2*x-1)/3.14-0.5)^K
	end
	# f(0) = 1.0 and f(1) = 0.
	return [(Tmax - Tmin)*f(t) + Tmin for t in range(0, stop = 1, length = nT)]
end

get_linear_cooling_schedule(Tmin, Tmax, nT) = return collect(reverse(range(Tmin, stop = Tmax, length = nT)))


