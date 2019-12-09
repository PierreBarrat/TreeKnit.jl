export compute_energy

let n::Int64=0
	global increment_n() = (n+=1)
	global reset_n() = (n=0)
	global get_n() = n
end

"""
"""
function compute_energy(conf::Array{Bool,1}, g::Graph)
	length(conf) != length(g.leaves) && error("`conf` and `g` do not have the same length")
	E = 0
	if length(g.leaves) + length(g.internals) == 1
		return E
	end
	for (i,s) in enumerate(conf)
		if s
			for k1 in 1:g.K
				# Ancestors in tree k1
				a1 = g.leaves[i].anc[k1]
				# If ancestor is identical to leaf for given configuration (i.e. only one spin up), go up
				# @time a1.conf[conf]
				while !a1.isroot && nspinup(a1.conf, conf) == 1
					a1 = a1.anc::SplitNode
				end
				for k2 in (k1+1):g.K
					# Same for 
					a2 = g.leaves[i].anc[k2]
					while !a2.isroot && nspinup(a2.conf, conf) == 1
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
					# if a1.conf[conf] and a2.conf[conf] are disjoint or if one contains the other, E=0
					# @time E += (a1.conf[conf] != a2.conf[conf] && !isnothing(findfirst(x->x==-1, a1.conf[conf] - a2.conf[conf])) &&  !isnothing(findfirst(x->x==-1, a2.conf[conf] - a1.conf[conf])) )
					if !are_disjoint(a1.conf, a2.conf, conf) && !is_contained(a1.conf, a2.conf, conf) && !is_contained(a2.conf, a1.conf, conf)
						E += 1
					end
				end
			end
		end
	end
	return E
end
function nspinup(nodeconf, conf)
	n = 0
	for i in 1:length(nodeconf)
		if nodeconf[i] && conf[i]
			n += 1
		end
	end
	return n
end
function are_disjoint(nconf1, nconf2, conf)
	for (i,s) in enumerate(conf)
		if s 
			if nconf1[i] === nconf2[i]
				return false
			end
		end
	end
	return true
end
function is_contained(nconf1, nconf2, conf) # is 1 in 2 ? 
	for (i,s) in enumerate(conf)
		if s 
			if nconf1[i] && !nconf2[i]
				return false
			end
		end
	end
	return true
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
	oconf = [copy(_conf)]
	for m in 1:M
		E, F = mcmcstep!(_conf, g, F, T, γ)
		ee[m+1] = E
		ff[m+1] = F
		# If new minimum is found
		if F < Fmin
			Fmin = F
			oconf = [copy(_conf)]
		end
		# If equal minimum is found
		if F == Fmin && mapreduce(x->x!=_conf, *, oconf)
			push!(oconf, copy(_conf))
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
	oconf = [ones(Bool, length(g.leaves))]
	E = [compute_energy(oconf[1],g)]
	F = Array{Float64,1}([E[1]])
	Fmin = F[1]
	for T in Trange
		tmp_oconf, e, f = SplitGraph.doMCMC(g, oconf[rand(1:length(oconf))], M, T=T,γ=γ)
		append!(E,e)
		append!(F,f)
		# If a better conf is found than all configurations in oconf (which is the min of `f` from doMCMC), completely replace oconf
		if findmin(f)[1] < Fmin
			oconf = tmp_oconf
			Fmin = findmin(F)[1]
		# If equally good confs have been found
		elseif findmin(f)[1] == Fmin
			append!(oconf, tmp_oconf)
			oconf = unique(oconf)
		end
	end
	return oconf,E,F
end

"""
"""
function count_mismatches(g::Graph)
	conf = ones(Bool, length(g.leaves))
	return compute_energy(conf, g)
end

"""
"""
function count_mismatches(t::Vararg{Tree})
	treelist = collect(t)
	mcc = maximal_coherent_clades(treelist)
	mcc_names = name_mcc_clades!(treelist, mcc)
	for (i,t) in enumerate(treelist)
		treelist[i] = reduce_to_mcc(t, mcc)
	end
	g = trees2graph(treelist)
	conf = ones(Bool, length(g.leaves))
	return compute_energy(conf, g)
end


