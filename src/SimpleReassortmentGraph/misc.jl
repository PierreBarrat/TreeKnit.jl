const rand_str_length = 8

make_random_label() = "ARGNode_$(randstring(rand_str_length))"
make_random_label(str::String) = str * "_" * randstring(rand_str_length)
make_random_label(i::Int) = "ARGNode_$(i)"
make_random_label(str, i) = str * "_" * i


function check(n::ARGNode; strict=true)
	#General topology
	## No ancestor ==> tau is missing
	@assert mapreduce(
		(a,τ) -> !isnothing(a) || (isnothing(a) && ismissing(τ)),
		&,
		n.anc, n.tau
	) "No ancestor and non-missing `tau` field (node $(n.label))"

	# Colors
	## color[i] <==> (root[i] or (ancestor[i] and the anc has color `i` too))
	## (i.e. not color --> root is false and ancestor is nothing)
	for (i,c) in enumerate(color(n))
		if c
			@assert isroot(n, i) || !isnothing(ancestor(n, i)) && hascolor(ancestor(n,i), i) "Inconsistent color / root / ancestor for $(label(n))"
		end
	end

	# If it's a leaf, it should have all colors
	if isleaf(n)
		@assert degree(n) == length(color(n))
	end

end

function Base.show(io::IO, n::ARGNode)
    if !get(io, :compact, false)
        nodeinfo(io, n)
    end
end
Base.show(n::ARGNode) = show(stdout, n)

"""
    nodeinfo(io, node)

Print information about `node`.
"""
function nodeinfo(io, node)
    print(io, "$(label(node)) ")
  	if node.hybrid
  		print(io, "(hybrid) ")
  	end
  	println(io, ":")
  	println(io, "Color: $(color(node))")
  	println(io, "Ancestors: $([label(x) for x in ancestors(node)])")
  	println(io, "Branch length: $(tau(node))")
  	println(io, "Children: $([label(x) for x in children(node)])")
end
