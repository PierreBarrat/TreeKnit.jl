"""
	write_mccs(outfile::Function, MCCs::AbstractDict[, segments])
	write_mccs(
		outfile::AbstractString,
		MCCs::AbstractDict,
		segments = get_segments(MCCs),
	)
	write_mccs(of, MCCs::AbstractArray)

Write MCCs to file.
"""
function write_mccs(
	outfile::AbstractString,
	MCCs::AbstractDict,
	segments = get_segments(MCCs),
)
	write_mccs((s1,s2) -> outfile, MCCs, segments)
end
function write_mccs(outfile::Function, MCCs::AbstractDict, segments = get_segments(MCCs))
	for i in 1:length(segments), j in (i+1):length(segments)
		write_mccs(outfile(segments[i], segments[j]), MCCs[segments[i], segments[j]])
	end

	return nothing
end

function write_mccs(of, MCCs::AbstractArray)
	open(of, "w") do w
		for m in MCCs
			for x in m[1:end-1]
				write(w, x * ",")
			end
			write(w, m[end] * "\n")
		end
	end

	return nothing
end

function get_segments(MCCs)
	segments = []
	for k in keys(MCCs)
		push!(segments, k...)
	end

	return unique(sort(segments))
end

function read_mccs(file)
	map(eachline(file)) do m
       String.(split(m, ','))
    end
end
