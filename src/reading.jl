"""
Read mutations from treetime built nexus file. Return a dictionnary "node label" => (mutations, ...)
"""
function parse_nexus(infile)
	f = readlines(infile)
	l = f[findfirst(x->x=="Begin Trees;", f) + 1]
	inodes = findall(x->x==':', l)
	out = Dict{Any,Any}()
	for i in inodes
		# Getting node name
		j = i
		while l[j]!=',' && l[j]!='(' && l[j]!=')'
			j -= 1
		end
		name = l[j+1:i-1]
		muts = Array{String,1}(undef, 0)
		# If there are mutations for this node, I will find a '[' before the next ',' or ')'
		j = i
		while l[j]!=',' && l[j] != ')' && l[j] != '(' && l[j] != ';' && l[j]!='['
			j += 1
		end
		# Getting mutations
		if l[j] == '['
			j2 = j
			while l[j2]!=']'
				j2 += 1
			end
			sub = l[j+1:j2-1]
			muts = String.(split(split(sub, '"')[2], ','))
		end
		out[name] = parse_muts(muts)
	end
	return out
end

