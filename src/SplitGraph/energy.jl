export compute_energy

"""
"""
function compute_energy(conf::Array{Bool,1}, g::Graph)
	length(conf) != length(g.leaves) && error("`conf` and `g` do not have the same length")
	E = 0
	for (i,s) in enumerate(conf)
		if s
			for k1 in 1:g.K
				# Ancestors in tree k1
				a1 = g.leaves[i].anc[k1]
				# If ancestor is identical to leaf for given configuration (i.e. only one spin up), go up
				# @time a1.conf[conf]
				while !a1.isroot && sum(a1.conf[conf]) == 1
					a1 = a1.anc
				end
				for k2 in (k1+1):g.K
					# Same for 
					a2 = g.leaves[i].anc[k2]
					while !a2.isroot && a2.conf[conf] == g.leaves[i].conf[conf]
						a2 = a2.anc
					end
					# Mismatch
					# println(a1.conf[conf] != a2.conf[conf] , " ", !isnothing(findfirst(x->x==-1, a1.conf[conf] - a2.conf[conf])) ||  !isnothing(findfirst(x->x==-1, a2.conf[conf] - a1.conf[conf])))
					# E += (a1.conf[conf] != a2.conf[conf])
					# if a1.conf[conf] != a2.conf[conf] 
					# 	println(a1.conf[conf])
					# 	println(a2.conf[conf])
					# 	println()
					# end
					E += (a1.conf[conf] != a2.conf[conf] && !isnothing(findfirst(x->x==-1, a1.conf[conf] - a2.conf[conf])) &&  !isnothing(findfirst(x->x==-1, a2.conf[conf] - a1.conf[conf])) )
				end
			end
		end
	end
	return E
end

"""
"""
function compute_F(conf::Array{Bool,1}, g::Graph, γ::Real) #this is a mu O_o
	E = compute_energy(conf, g)
	return E + γ*(length(conf) - sum(conf))
end

"""
"""
function doMCMC(g::Graph, conf::Array{Bool,1}, M::Int64; T=1, γ=1)
	_conf = copy(conf)
	E = compute_energy(_conf, g)
	F = E + γ*(length(conf) - sum(conf))
	## 
	ee = zeros(Real,M+1)
	ff = zeros(Real,M+1)
	ee[1] = E
	ff[1] = F
	## 
	Fmin = F
	oconf = copy(_conf)
	for m in 1:M
		E, F = mcmcstep!(_conf, g, F, T, γ)
		ee[m+1] = E
		ff[m+1] = F
		if F < Fmin
			Fmin = F
			oconf = copy(_conf)
		end
	end
	return oconf,ee,ff
end

"""
"""
function mcmcstep!(conf, g, F, T, γ)
	i = rand(1:length(conf))
	conf[i] = !conf[i]
	Enew = compute_energy(conf, g)
	Fnew = Enew + γ*(length(conf) - sum(conf))
	if Fnew < F || exp(-(Fnew-F)/T) > rand()
		return Enew, Fnew
	else
		conf[i] = !conf[i]
		return (round(Int64, F-γ*(length(conf) - sum(conf))), F)
	end
end

"""
"""
function sa_opt(g::Graph ; Trange=1.:-0.01:0.1, γ=1.05, M=1000)
	confinit = ones(Bool, length(g.leaves))
	E = []
	F = []
	for T in Trange
		confinit, e, f = SplitGraph.doMCMC(g, confinit, M, T=T,γ=γ)
		append!(E,e)
		append!(F,f)
	end
	return confinit,E,F
end