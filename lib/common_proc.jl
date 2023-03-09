module common_sub
using CSV
using DataFrames
using DataFramesMeta

export get_short_area_code

function process_area_name(list_area_name)
    tbr1 = Vector{String}(undef, length(list_area_name))
    tbr2 = Vector{String}(undef, length(list_area_name))
  
    for i in eachindex(list_area_name)
      area_name = list_area_name[i]
      ch = findfirst(",", area_name)
      ch2 = findfirst("-", area_name)
      if ch2 === nothing
        tbr2[i] = area_name
      else
        tbr2[i] = first(area_name, ch2[1] - 1)
  
      end
      if ch === nothing
        tbr1[i] = area_name
      else
        tbr1[i] = first(area_name, ch[1] - 1)
      end
  
    end
    return ([tbr1 tbr2])#, new_name_2)
  end
  
  
  function get_short_area_code(df_names::DataFrame)

    transform!(df_names, :, :area_name => process_area_name => [:short1, :brief_area])
    return (df_names[:, [:area_code, :brief_area]])
  
  end
  
  
end