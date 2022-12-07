"""
	write_mccs(of, MCCs::AbstractArray)

Write MCCs to file.
"""
function write_mccs(of, MCCs::AbstractArray, mode="w")
	open(of, mode) do w
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

