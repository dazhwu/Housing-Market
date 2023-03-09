
__precompile__(true)

module vector_utilities

export get_gaps, consecutive_ranges

function get_gaps(vec::Vector{Int})
    len=length(vec)
    lowerBounds=[vec[i]+1 for i =1:len-1]
    upperBounds=[vec[i]-1 for i=2:len]
    mask=lowerBounds .<= upperBounds
    return (lowerBounds[mask], upperBounds[mask])
end

function consecutive_ranges(vec::Vector{Int})
    len = length(vec)
    groups = Vector{Vector{eltype(vec)}}()
  
    i = j = 1
    while i <= len && j <= len
      j = i
      while j < len  && vec[j] + 1 == vec[j + 1] 
        j += 1
      end
      push!(groups,vec[i:j])
      i = j + 1
    end
  
    return groups
  end

end
# lb, ub=get_gaps([3,9,10,11,12,26,30,38,39])
# print(lb)
# print(ub)
# print(consecutive_ranges([3,9,10,11,12,26,30,38,39]))