"""
	branch_likelihood(n1::Real, n2::Real, μ1::Real, μ2::Real)

Given that `n1` and `n2` mutations were observed on two branches, with mutation rates `μ1` and `μ2`, return the log-ratio of:
- the probability that the two branch lengths `t1` and `t2` are the same and equal to `(n1 + n2)/(μ1 + μ2)`.
- the probability that the two branch lenghts are different and equal to respectively `n1/μ1` and `n2/μ2`. 
"""
function branch_likelihood(n1::Real, n2::Real, μ1::Real, μ2::Real)
	L = 0.
	n1 != 0. && (L += n1*log(μ1/(μ1+μ2) * (n1+n2)/n1))
	n2 != 0. && (L += n2*log(μ2/(μ1+μ2) * (n1+n2)/n2))
	return L
end



function conf_likelihood(conf::Array{Bool,1}, g::Graph, μ, trees; v=false, mode=:mutations)
	if in(mode, (:mutations, :nmut))
		return conf_likelihood_(TreeTools.node_diffmut, conf, g, μ, trees, v=v)
	elseif in(mode, (:time, :times, :divtime))
		return conf_likelihood_(TreeTools.node_divtime, conf, g, μ, trees, v=v)
	else
		@error "Unrecognized mode $(mode)."
	end
end

"""
	conf_likelihood_(divfunction, conf::Array{Bool,1}, g::Graph, μ, trees; v=false)

For each leaf `n` remaining in `conf`, find the first non-trivial ancestors in graph `g`. For each pair of ancestors `a1` and `a2`, compare the likelihood  that the real branch length from `n` to `a1` and `a2` is the same, with the likelihood that it is different. 
"""
function conf_likelihood_(divfunction, conf::Array{Bool,1}, g::Graph, μ, trees; v=false)
	L = 0.
	if length(g.leaves) + length(g.internals) == 1
		return L
	end
	for (i,s) in enumerate(conf)
		if s
			# Similar to compute_energy 
			for k1 in 1:g.K
				# Tree node of k1 corresponding to graph leaf i
				tn1 = trees[k1].lleaves[g.labels[i]]
				# Going up in SplitGraph and tree at the same time
				a1 = g.leaves[i].anc[k1]
				ta1 = tn1.anc
				while onespinup(a1.conf, conf) && !a1.isroot 
					a1 = a1.anc::SplitNode
					ta1 = ta1.anc
				end
				for k2 in (k1+1):g.K
					# Tree node of k2 corresponding to graph leaf i
					tn2 = trees[k2].lleaves[g.labels[i]]
					# Going up in SplitGraph and tree at the same time
					a2 = g.leaves[i].anc[k2]
					ta2 = tn2.anc
					while onespinup(a2.conf, conf) && !a2.isroot 
						a2 = a2.anc
						ta2 = ta2.anc
					end
					# If a1.conf and a2.conf are equal, we just inferred that this branch is common to trees k1 and k2
					# if are_equal(a1.conf, a2.conf, conf)
					if are_equal_with_resolution(g, a1.conf, a2.conf, conf, k1, k2)
						n1 = divfunction(tn1, ta1)
						n2 = divfunction(tn2, ta2)
						v && println("No inconsistency for leaf $i ($(g.labels[i]))")
						v && println("Inferring that times $(n1/μ[k1]) on tree $k1 and $(n2/μ[k2]) on tree $k2 are the same.")
						v && println("Corresponding number of mutations: $n1 with rate $(μ[k1]) and $n2 with rate $(μ[k2]).")
						dL = branch_likelihood(n1, n2, μ[k1], μ[k2])

						# in dL, μ[k1]*τ1 should ultimately be replaced by n1 which is the observed number of mutations (same goes for 2)
						# What is currently there only makes sense for artificial data, where t is known exactly. 
						L += dL
						v && println("Corresponding likelihood is $dL.")
					end
				end
			end
		end
	end
	return L
end
"""
"""
function conf_likelihood_times(conf::Array{Bool,1}, g::Graph, μ, trees; v=false)
	L = 0.
	for (i,s) in enumerate(conf)
		if s
			# Similar to compute energy
			for k1 in 1:g.K
				# Tree node of k1 corresponding to graph leaf i
				tn1 = trees[k1].lleaves[g.labels[i]]
				# Going up in SplitGraph and tree at the same time
				a1 = g.leaves[i].anc[k1]
				ta1 = tn1.anc
				while onespinup(a1.conf, conf) && !a1.isroot 
					a1 = a1.anc::SplitNode
					ta1 = ta1.anc
				end
				for k2 in (k1+1):g.K
					# Tree node of k2 corresponding to graph leaf i
					tn2 = trees[k2].lleaves[g.labels[i]]
					# Going up in SplitGraph and tree at the same time
					a2 = g.leaves[i].anc[k2]
					ta2 = tn2.anc
					while onespinup(a2.conf, conf) && !a2.isroot 
						a2 = a2.anc
						ta2 = ta2.anc
					end
					# If a1.conf and a2.conf are equal, we just inferred that this branch is common to trees k1 and k2
					if are_equal(a1.conf, a2.conf, conf)
						τ1 = TreeTools.node_divtime(tn1, ta1)
						τ2 = TreeTools.node_divtime(tn2, ta2)
						dL = branch_likelihood(μ[k1]*τ1, μ[k2]*τ2, μ[k1], μ[k2])
						# in dL, μ[k1]*τ1 should ultimately be replaced by n1 which is the observed number of mutations (same goes for 2)
						# What is currently there only makes sense for artificial data, where t is known exactly. 
						L += dL
						v && println("No inconsistency for leaf $i ($(g.labels[i]))")
						v && println("Inferring that times $(τ1) on tree $k1 and $(τ2) on tree $k2 are the same.")
						v && println("Corresponding likelihood is $dL.")
					end
				end
			end
		end
	end
	return L
end
