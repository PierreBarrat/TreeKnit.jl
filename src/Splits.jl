export Split

struct Split{T}
	dat::Tuple{Set{T}, Set{T}}
end
Split(x::Set{T}, y::Set{T}) where T = return Split((x,y))
Split(x::T, y::Set{T}) where T = return Split((Set([x]), y))
Split(x::Set{T}, y::T) where T = return Split((x, Set([y])))
Split(x::T, y::T) where T = return Split((Set([x]), Set([y])))
Split(x::T) where T = return Split(Set([x]), Set{T}())
Split(x::Set{T}) where T = return Split(x, Set{T}())

getindex(s::Split{T}, i::Int64) where T = s.dat[i]
setindex!(s::Split{T}, x::Set{T}, i::Int64) where T = (s.dat[i] = x)
function ==(s::Split{T}, t::Split{T}) where T 
	if !isempty(s[1]) && (s[1] == t[1] || s[1] == t[2])
		return true
	elseif !isempty(s[2]) && (s[2] == t[1] || s[2] == t[2])
		return true
	else
		return false
	end
end


function Split(n::TreeNode)
	return Split(Set(TreeTools.node_leavesclade_labels(n)))
end