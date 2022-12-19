"""
	write(io::IO, arg::ARG)
	write(filename::AbstractString, arg::ARG)

Write `arg` using the extended Newick format.
"""
write(filename::AbstractString, arg::ARG) = write(filename, extended_newick(arg))
write(io::IO, arg::ARG) = write(io, extended_newick(arg))

#=
Strategy
(i) If the global root is shared
Go recursively down this root and write a newick.
If a hybrid node is met for the first time, store it in a list of hybrid nodes met
If it is met for the second time, stop recursing down
(ii) none of the roots are shared
feed the function with each root in turn (and the list of hybrids)
(iii) One of the root is shared, not the other one
start with the non shared one. should work allright
=#
function extended_newick(arg::ARG)
	hybrids = Dict()

	if arg.roots[1] == arg.roots[2]
		# (i)
		str = extended_newick!(arg.roots[1], nothing, hybrids)
		return str * ";"
	else
		if !isshared(arg.roots[1]) && !isshared(arg.roots[2])
			#(ii)
			str1 = extended_newick!(arg.roots[1], nothing, hybrids)
			str2 = extended_newick!(arg.roots[2], nothing, hybrids)
			return "($(str1),$(str2))GlobalRoot[&segments={0,1}]:0.;"
		elseif isshared(arg.roots[1])
			#(iii)
			str = extended_newick!(arg.roots[2], nothing, hybrids)
			return str * ";"
		elseif isshared(arg.roots[2])
			str = extended_newick!(arg.roots[1], nothing, hybrids)
			return str * ";"
		end
	end
end

function extended_newick!(a_r::ARGNode, anc, hybrids::Dict)
	# println("$(label(a_r))")
	nwk = ""
	# Dealing with children of a_r
	if !isleaf(a_r) && !haskey(hybrids, label(a_r))
		nwk *= "("
		for c in children(a_r)
			nwk *= extended_newick!(c, a_r, hybrids)
			nwk *= ","
		end
		nwk = nwk[1:end-1] # Removing trailing ','
		nwk *= ")"
	end
	# Adding a_r to hybrid if necessary
	if ishybrid(a_r) && !haskey(hybrids, label(a_r))
		i = length(hybrids) +1
		hybrids[label(a_r)] = "#H$(i)"
	end
	# Writing a_r itself
	a_r_label = if ishybrid(a_r)
		label(a_r) * hybrids[label(a_r)]
	else
		label(a_r)
	end
	nwk *= a_r_label
	nwk *= nwk_node_data(a_r, anc)

	return nwk
end

"""
	nwk_node_data(n::ARGNode, a)

Return a string of the form "[&segments=0,1]:$(time)".

Different cases to choose from:
	(i) Node is shared, not a hybrid, not the root of only one tree
  its two times should be the same
  [&segments={0,1}]
	(ii) Node is not shared
  It should only have one color
  [&segments={c}]
	(iii) Node is shared and a hybrid
  We need to know what the color of the branch we arrived from is
  [&segments={c}]
	(iv) Node is the root of one tree and not of the other
  If `c` is the non-root color
  [&segments={c}]
"""
function nwk_node_data(n::ARGNode, a)
	if isshared(n) && !ishybrid(n) && !is_partial_root(n)
		# (i)
		τ1 = branch_length(n)[1]
		τ2 = branch_length(n)[2]
		@assert (ismissing(τ1) && ismissing(τ2)) || isapprox(τ1, τ2)
		return nwk_node_data(τ1, 1, 2)

	elseif isshared(n) && ishybrid(n)
		#(iii)
		ec = edgecolor(a, n)
		if ec[1] && !ec[2]
			return nwk_node_data(branch_length(n)[1], 1)
		elseif !ec[1] && ec[2]
			return nwk_node_data(branch_length(n)[2], 2)
		elseif ec[1] && ec[2]
			@assert isapprox(branch_length(n)[1], branch_length(n)[2])
			return nwk_node_data(branch_length(n)[1], 1, 2)
		end

	elseif isshared(n) && isroot(n)
		# (iv)
		for c in 1:2
			if !isroot(n, c)
				@assert isroot(n, othercolor(c))
				return nwk_node_data(branch_length(n)[c], c)
			end
		end
		@error "$n"
	elseif !isshared(n)
		# (ii)
		for c in 1:2
			if hascolor(n, c)
				@assert !hascolor(n, othercolor(c))
				return nwk_node_data(branch_length(n)[c], c)
			end
		end
		@error "No color for node $n"
	else
		@error "Unclear state for node $n"
	end
end


function nwk_node_data(τ, c::Int)
	ismissing(τ) ? "[&segments={$(c-1)}]" : "[&segments={$(c-1)}]:$(τ)"
end
function nwk_node_data(τ, c1::Int, c2::Int)
	clr = "$(min(c1,c2)-1),$(max(c1,c2)-1)"
	ismissing(τ) ? "[&segments={$(clr)}]" : "[&segments={$(clr)}]:$(τ)"
end
