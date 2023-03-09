using CSV
using DataFrames
using DataFramesMeta
using ShiftedArrays
using Statistics
using LinearAlgebra
using Distributions
using StateSpaceModels
using HypothesisTests
using DelimitedFiles
using Plots
using LaTeXStrings
using Downloads

include("./lib/vector_utilities.jl")
include("./lib/tvgc.jl")
include("./lib/ivx.jl")
include("./lib/exuber_old.jl")
include("./lib/common_proc.jl")
include("./lib/spillover/SpillPap.jl")

df_names = common_sub.get_short_area_code(CSV.read("cu_area.txt", DataFrame))

df0 = CSV.read("combined.csv", DataFrame)
df0 = innerjoin(df0, df_names, on=:area_code)
df0.strSeq = string.(df0.Year) .* "Q" .* string.(df0.Quarter)

strSeq = sort(unique(df0.strSeq))  #[3:end]


#df1=CSV.read("bubble.csv", DataFrame)
@transform!(df0, :y = :log_index_nsa - :log_rent)
df1 = df0[!, [:y, :brief_area, :seq, :cpi, :log_emp, :log_rent, :une_rate, :log_population, :m_rate, :log_Close]]  # :log_num_permits
dropmissing!(df1)
num_rows, _ = size(df1)
MSAs = unique(df1.brief_area)
#df_names[df_names.brief_area.==MSAs[i], :short2][1]
num_MSAs = length(MSAs)

average_length = Int(num_rows / num_MSAs) - 2
wide_table_F = Matrix{Float64}(undef, average_length, num_MSAs)
wide_table_B = Matrix{Float64}(undef, average_length, num_MSAs)
log_price_to_rent = Matrix{Float64}(undef, Int(num_rows / num_MSAs), num_MSAs)
seq = Vector{Float64}(undef, average_length)



for index = eachindex(MSAs)

  msa = MSAs[index]
  partial = df1[(df1.brief_area.==msa), [:y, :seq, :cpi, :log_emp, :une_rate, :log_rent, :log_population, :m_rate, :log_Close]]  #:log_num_permits

  num_obs, num_vars = size(partial)
  if index == 1
    seq = Vector(partial.seq)[2:num_obs]

  end
  temp_y = Vector(partial.y)  # 1 to num_obs  
  temp_x = Matrix(partial[:, 3:num_vars])
  log_price_to_rent[:, index] = temp_y
  AR = true
  wide_table_F[:, index], wide_table_B[:, index] = _ivx.decompose(temp_y, temp_x, AR, 4, 1, 0.3, 0.02)
end

df_y = DataFrame(wide_table_B, :auto)
rename!(df_y, Symbol.(MSAs))
insertcols!(df_y, 1, :strSeq => strSeq[3:end])

CSV.write("bubble.csv", df_y)

# df_y = DataFrame(log_price_to_rent, :auto)
# rename!(df_y, Symbol.(MSAs))
# insertcols!(df_y, 1, :strSeq => strSeq)

#T, col_width = size(log_price_to_rent)
T, col_width = size(wide_table_B)

p = 2
typ = "Const"

est =SpillPap.varEstimate(wide_table_B, p, typ)
t = SpillPap.spilloverTable(est, 8; fevd_function="genFEVD", nocorr = false)
tables=SpillPap.SpilloverTable(t)

SpillPap.printTable(tables, MSAs)

trend = false
IC = "bic"

max_lag = 4
window_size = floor(Int64, 0.2 * T)  #18
# 12              #% 12 years
robust = true
alpha = 0.95
control_size = 6  #12

results = _tvgc.tvgc(df_y, trend, max_lag, window_size, control_size, alpha, robust, IC)

num_plots = size(results[5], 1)

the_folder = "./png/TVGC/recursive/"
rm(the_folder, force=true, recursive=true)
mkdir(the_folder)

for i in 1:num_plots
  the_label = results[5][i, 1] * "\$\\longmapsto\$" * results[5][i, 2]
  plot(strSeq[window_size:T], results[1][:, i], label=the_label, dpi=600, lw=2, linecolor=:black)
  plot!(strSeq[window_size:T], results[2][:, i], label="95% cv", line=(:dot, 2))
  savefig(the_folder* results[5][i, 1] * " to " * results[5][i, 2] * ".png")
end

the_folder = "./png/TVGC/rolling window/"
rm(the_folder, force=true, recursive=true)
mkdir(the_folder)

for i in 1:num_plots
  the_label = results[5][i, 1] * "\$\\longmapsto\$" * results[5][i, 2]
  plot(strSeq[window_size:T], results[3][:, i], label=the_label, dpi=600, lw=2, linecolor=:black)
  plot!(strSeq[window_size:T], results[4][:, i], label="95% cv", line=(:dot, 2))
  savefig(the_folder* results[5][i, 1] * " to " * results[5][i, 2] * ".png")
end

