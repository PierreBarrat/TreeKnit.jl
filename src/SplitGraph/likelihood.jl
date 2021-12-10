"""
	branch_likelihood(t1, t2, L1, L2)

Given two branches of length `t1` and `t2`, with segments of length `L1` and `L2`, return the log-ratio of:
- the probability that the two branch lengths `t1` and `t2` are the same and equal to `(t1*L1 + t2*L2) / (L1 + L2)`.
- the probability that the two branch lenghts are different and equal to respectively `t1` and `t2`. 
"""
function branch_likelihood(t1::Real, t2::Real, L1, L2) 
	n1s = (t1*L1 + t2*L2) / (L1 + L2) * L1
	n2s = (t1*L1 + t2*L2) / (L1 + L2) * L2
	n1 = t1*L1
	n2 = t2*L2
	L = 0.
	if n1 != 0 
		L += n1 - n1s + n1 * log(n1s/n1)
	else
		L += -n1s
	end
	if n2 != 0 
		L += n2 - n2s + n2 * log(n2s/n2)
	else
		L += -n2s
	end
	return L 
end

branch_likelihood(t1::Missing, t2::Real, L1, L2) = missing
branch_likelihood(t1::Real, t2::Missing, L1, L2) = missing
branch_likelihood(t1::Missing, t2::Missing, L1, L2) = 0.
branch_likelihood(t1, t2, L1, L2) = 0.


function conf_likelihood(conf::Array{Bool,1}, g::Graph, μ, trees; v=false, mode=:mutations)
	conf_likelihood_(TreeTools.divtime, conf, g, μ, trees, v=v)
end

"""
	conf_likelihood_(divfunction, conf::Array{Bool,1}, g::Graph, μ, trees; v=false)

For each leaf `n` remaining in `conf`, find the first non-trivial ancestors in graph `g`. For each pair of ancestors `a1` and `a2`, compare the likelihood  that the real branch length from `n` to `a1` and `a2` is the same, with the likelihood that it is different. 
"""
function conf_likelihood_(divfunction, conf::Array{Bool,1}, g::Graph, seq_lengths, trees; v=false)
	L = 0.
	Z = 0.
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
				while is_trivial_split(a1.conf, conf) && !a1.isroot
					a1 = a1.anc::SplitNode
					ta1 = ta1.anc
				end
				for k2 in (k1+1):g.K
					# Tree node of k2 corresponding to graph leaf i
					tn2 = trees[k2].lleaves[g.labels[i]]
					# Going up in SplitGraph and tree at the same time
					a2 = g.leaves[i].anc[k2]
					ta2 = tn2.anc
					while is_trivial_split(a2.conf, conf) && !a2.isroot
						a2 = a2.anc
						ta2 = ta2.anc
					end
					# If a1.conf and a2.conf are equal, we just inferred that this branch is common to trees k1 and k2
					if (!get_resolve() && are_equal(a1.conf, a2.conf, conf)) ||
						(get_resolve() && are_equal_with_resolution(g, a1.conf, a2.conf, conf, k1, k2))
						tau1 = divfunction(tn1, ta1)
						tau2 = divfunction(tn2, ta2)
						# v && println("No inconsistency for leaf $i ($(g.labels[i]))")
						# v && println("Inferring that times $(tau1) on tree $k1 and $(tau2) on tree $k2 are the same.")
						dL = branch_likelihood(tau1, tau2, seq_lengths[k1], seq_lengths[k2])
						L += dL
						Z += 1.
						# v && println("Corresponding likelihood is $dL.")
					end
				end
			end
		else # branch above `i` is not shared, since we removes `i` from the trees
			for k1 in 1:g.K
				# Tree node of k1 corresponding to graph leaf i
				tn1 = trees[k1].lleaves[g.labels[i]]
				ta1 = tn1.anc
				for k2 in 1:g.K
					# Tree node of k2 corresponding to graph leaf i
					tn2 = trees[k2].lleaves[g.labels[i]]
					ta2 = tn2.anc
					# Removing log-likelihood ratio, since it's L(shared) / L(non shared)
					tau1 = divfunction(tn1, ta1)
					tau2 = divfunction(tn2, ta2)
					dL = -branch_likelihood(tau1, tau2, seq_lengths[k1], seq_lengths[k2])
					L += dL
					Z += 1.
				end
			end
		end
	end
	# v && println("--> Final likelihood $(L/max(Z,1))")
	return L/max(Z,1)
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
				while is_trivial_split(a1.conf, conf) && !a1.isroot
					a1 = a1.anc::SplitNode
					ta1 = ta1.anc
				end
				for k2 in (k1+1):g.K
					# Tree node of k2 corresponding to graph leaf i
					tn2 = trees[k2].lleaves[g.labels[i]]
					# Going up in SplitGraph and tree at the same time
					a2 = g.leaves[i].anc[k2]
					ta2 = tn2.anc
					while is_trivial_split(a2.conf, conf) && !a2.isroot
						a2 = a2.anc
						ta2 = ta2.anc
					end
					# If a1.conf and a2.conf are equal, we just inferred that this branch is common to trees k1 and k2
					if are_equal(a1.conf, a2.conf, conf)
						τ1 = TreeTools.divtime(tn1, ta1)
						τ2 = TreeTools.divtime(tn2, ta2)
						dL = branch_likelihood(μ[k1]*τ1, μ[k2]*τ2, μ[k1], μ[k2])
						# in dL, μ[k1]*τ1 should ultimately be replaced by n1 which is the observed number of mutations (same goes for 2)
						# What is currently there only makes sense for artificial data, where t is known exactly. 
						L += dL
						v && @info "No inconsistency for leaf $i ($(g.labels[i]))"
						v && @info "Inferring that times $(τ1) on tree $k1 and $(τ2) on tree $k2 are the same."
						v && @info "Corresponding likelihood is $dL."
					end
				end
			end
		end
	end
	return L
end
