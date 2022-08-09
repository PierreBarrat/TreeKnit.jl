"""
	write_mccs(of, MCCs::AbstractArray)

Write MCCs to file.
"""
function write_mccs(filePath, MCCs::AbstractArray, mode="w")
	file_type = split(filePath, ".")[2]
	if file_type !="dat"
		print("File type not supported")
	end
	open(filePath, mode) do w
		for (i,m) in enumerate(MCCs)
			for x in m[1:end-1]
				write(w, x * ",")
			end
			if i < length(MCCs)
				write(w, m[end] * "\n")
			else
				write(w, m[end])
			end
		end
	end

	return nothing
end


function read_mccs(file)
	map(eachline(file)) do m
       String.(split(m, ','))
    end
end


function write_mccs(filePath, MCCs::MCC_set, mode="w")
	file_type = split(filePath, ".")[2]
	if file_type !="json"
		print("File type not supported")
	end

	open(filePath, mode) do w
		write(w, "{ \"MCC_dict\" : {\n")
		tree_pairs, Ms = iter_pairs(MCCs)
		for i in 1:binomial(MCCs.no_trees, 2)
			write(w, "\""*string(i)*"\": {\n \"trees\":[")
			s = sort(tree_pairs[i], lt=clt)
			write(w, "\""*string(s[1])*"\", \""*string(s[2])*"\"],")
			write(w, "\n\"mccs\": [")
			m = Ms[i]
			for x in m[1:end-1]
				write(w, string(x) * ",\n")
			end
			if i < length(Ms)
				write(w, string(m[end]) *"]\n},\n")
			else
				write(w, string(m[end])*"]\n}\n}\n}")
			end
		end
	end

	return nothing
end

