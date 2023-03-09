using CSV
using DataFrames
using Dates
using MarketData
using DataFramesMeta
#using DuckDB
using Statistics
using TimeSeries
using Plots
using ARCHModels
#using Downloads
using StringDistances

include("./lib/vector_utilities.jl")
include("./lib/tvgc.jl")
include("./lib/ivx.jl")
include("./lib/exuber_old.jl")
include("./lib/exuber.jl")
include("./lib/common_proc.jl")
include("./lib/spillover/SpillPap.jl")


function match_names(df1, df2, key1, key2)
	matched = String[]
	for s in df1[:, key1]
		_, ind = findnearest(s, df2[:, key2], Levenshtein())
		push!(matched, df2[ind, :area_code])
	end

	matched = DataFrame(area_code = matched, SP_area = df1[:, :SP_area])
	area = innerjoin(df2, matched, on = :area_code)
	return (area)
end

function download_employment(area_unemp::DataFrame)

	area = @chain CSV.read(download("https://download.bls.gov/pub/time.series/sm/sm.area"), DataFrame; normalizenames = true) begin
		@rsubset contains(:area_name, r"[,]")
		#@rsubset contains(:area_name, r"," )        
		@rtransform :brief = SubString(:area_name, 1, findfirst(",", :area_name)[1] - 1)
		@rtransform :area_code = string(:area_code)

	end

	area = match_names(area_unemp, area, :area_text, :brief)

	Emp = @chain CSV.read(download("https://download.bls.gov/pub/time.series/sm/sm.data.1.AllData"), DataFrame; normalizenames = true) begin

		@rtransform :area_code = SubString(:series_id, 6, 10)
		@rtransform :data_type = SubString(:series_id, 19, 20)
		@rsubset :data_type == "01"
		@rtransform :month = parse(Int32, SubString(:period, 2, 3))
		@rsubset :month <= 12
		@rtransform :industry_code = SubString(:series_id, 11, 18)
		@rsubset :industry_code == "00000000"
		@rtransform :season = SubString(:series_id, 3, 3)
		@rsubset :season == "S"
		@rtransform :date = Date(:year, :month, 1)
		@rtransform :value = parse(Float64, :value)
		@select(:date, :value, :area_code, :data_type)
	end

	Emp = innerjoin(Emp, area, on = "area_code")
	tbr = unstack(Emp, [:date, :SP_area], :data_type, :value)
	DataFramesMeta.rename!(tbr, [:date, :SP_area, :emp])
	return (tbr)
end

function download_unemployment(df_MSA)
	#unemployment data include state level，MSA level， county level
	area = @chain CSV.read(download("https://download.bls.gov/pub/time.series/la/la.area"), DataFrame; normalizenames = true) begin
		@rsubset :area_type_code == "B"
		@rtransform :brief = SubString(:area_text, 1, findfirst(",", :area_text)[1] - 1)

		@rtransform :brief = !isnothing(findfirst("-", :brief)) && SubString(:brief, 1, findfirst("-", :brief)[1] - 1)
		@select(:area_code, :area_text, :brief)
	end

	area = match_names(df_MSA, area, :SP_area, :brief)


	data_files = ["90-94", "95-99", "00-04", "05-09", "10-14", "15-19", "20-24"]
	combined = DataFrame(area_code = String[], measure_code = String[], date = Date[], value = Float64[])

	for f in eachindex(data_files)
		source = "https://download.bls.gov/pub/time.series/la/la.data.0.CurrentU" * data_files[f]
        #temp_df=CSV.read(download(source), DataFrame; normalizenames = true, types=[String, Int32, String, String, String])
		#CSV.write( data_files[f]*".csv" , temp_df)
        #println(size(temp_df))
		df = @chain CSV.read(download(source), DataFrame; normalizenames = true, types=[String, Int32, String, String, String]) begin
			@rtransform :area_code = SubString(:series_id, 4, 18)
			@rtransform :measure_code = SubString(:series_id, 19, 20)
			@rtransform :month = parse(Int32, SubString(:period, 2, 3))
			@rsubset :month <= 12		
			@rsubset contains(:value, r"-?[0-9]\d*(\.\d+)?$")	
			@rtransform :value=parse(Float64, :value)
			#@rsubset isa(:value, Number)
			@rtransform :date = Date(:year, :month, 1)
			@select(:area_code, :measure_code, :date, :value)
			
			#@rtransform :n=isa(:value, Number)
			#@select(:area_code, :measure_code, :date, :value)

		end
		
        println(data_files[f])
        println(size(df))
		append!(combined, df)
	end
	combined = innerjoin(combined, area, on = :area_code)
	combined = combined[:, [:SP_area, :date, :value, :measure_code]]

	temp_df=@rsubset(combined, :measure_code=="03")
	unemp_rate=unstack(temp_df, :date, :SP_area, :value)

	temp_df=@rsubset(combined, :measure_code=="04")
	unemp=unstack(temp_df, :date, :SP_area, :value)

	temp_df=@rsubset(combined, :measure_code=="05")
	emp=unstack(temp_df, :date, :SP_area, :value)

	temp_df=@rsubset(combined, :measure_code=="06")
	labor_force=unstack(temp_df, :date, :SP_area, :value)

	#measure_code	measure_text
	# 03	unemployment rate
	# 04	unemployment
	# 05	employment
	# 06	labor force

	#DataFramesMeta.rename!(tbr, [:date, :SP_area, :unemployment_rate, :unemployment, :employment, :labor_force])
	return (unemp_rate, unemp, emp, labor_force, area)

end

function download_price_index(MSA_CU::DataFrame)
	
	
	Price = @chain CSV.read(download("https://download.bls.gov/pub/time.series/cu/cu.data.0.Current"), DataFrame; normalizenames = true) begin

		@transform(:area_code = SubString.(:series_id, 5, 8))
		@rtransform(:item_code = SubString(:series_id, 9, length(:series_id)))
		@rtransform(:seasonal = SubString(:series_id, 3, 3))
		@rtransform :freq = SubString(:series_id, 4, 4)
		@rtransform :month = parse(Int32, SubString(:period, 2, 3))
		@rsubset contains(:item_code, r"(SEHA)|(SA0L2)")
		@rsubset :freq == "R"
		@rsubset :month <= 12
		@rsubset contains(:area_code, r"S[1-4]00")		
		@rtransform :date = Date(:year, :month, 1)
		#@rsubset contains(:area_code, r"S[1-4][0-9][0A-Z]")
		@select(:area_code, :date, :item_code, :value)
	end

	Price=innerjoin(Price, MSA_CU, on=:area_code)

	temp_rent = @chain Price begin
		@rsubset contains(:item_code, r"SEHA")
	end

	rent=unstack(temp_rent, :date, :SP_area, :value)

	temp_cpi = @chain Price begin
		@rsubset contains(:item_code, r"SA0L2")
	end
	cpi=unstack(temp_cpi, :date, :SP_area, :value)


	return (rent, cpi)
end

function download_SP_housing_index()

	#"DA", Dallas
	list_prefix = ["SF", "LX", "SE", "SD", "PH", "NY", "CH", "MI", "BO", "DN", "WD", "AT", "TP", "LV", "PO", "CR", "MN", "DE", "CE"]
	list_name = ["San Francisco", "Los Angeles", "Seattle", "San Diego", "Phoenix", "New York", "Chicago", "Miami", "Boston", "Denver", "Washington", "Atlanta", "Tampa", "Las Vegas", "Portland", "Charlotte", "Minneapolis", "Detroit", "Cleveland"]

	list_postfix = ["XRSA", "XRNSA"]
	list_postname = ["SPCase SA", "SPCase NonSA"]

	df = DataFrame(date = Date[], SP_area = String[], Index_type = String[], value = Float64[])
	for j in 1:length(list_postfix)
		for i in 1:length(list_prefix)
			temp = fred(list_prefix[i] * list_postfix[j])
			println(list_prefix[i] * list_postfix[j])
			temp_df = DataFrame(date = timestamp(temp), SP_area = list_name[i], Index_type = list_postfix[j], value = values(temp))
			append!(df, temp_df)
		end
	end
	return (df)
end

function download_mortgage_rate()
	mortgage = DataFrame(fred("MORTGAGE30US"))
	monthly_rate = @chain mortgage begin
		@transform(:month = Dates.month.(:timestamp))
		@transform(:year = Dates.year.(:timestamp))
		@by([:year, :month], :mortgage = mean(:VALUE))
        @rtransform :date=Date(:year, :month, 1)
        @select(:date, :mortgage)
	end
	return (monthly_rate)
end



function DY12(wide_var::DataFrame)
	n_col = size(wide_var, 2)
	data::Matrix{Float64} = Matrix(wide_var[:, 2:n_col])
	est = SpillPap.varEstimate(data, 1, "Const")
	t = SpillPap.spilloverTable(est, 8; fevd_function = "genFEVD", nocorr = false)
	tables = SpillPap.SpilloverTable(t)
	return (tables)
end

function diff(wide_table::DataFrame, order::Int64=1, ln::Bool=false)
	
	mat= Matrix(wide_table[:, 2:ncol(wide_table)])
	diff_mat = Matrix{Float64}(undef, nrow(wide_table)-order, ncol(wide_table)-1)

	if ln
		mat=log.(mat)
	end

	diff_mat=	@views mat[1+order:end, :] - mat[1:end-order, :]
	

	list_date = Vector(wide_table[:, 1])[1+order:end]

	df = DataFrame(diff_mat, :auto)

	DataFramesMeta.rename!(df, Symbol.(names(wide_table)[2:end]))
	insertcols!(df, 1, :date => list_date)

	return (df)
end

function format_date(d::Date, freq::String)
	y=Dates.year(d)
	m=Dates.month(d)
	d=Dates.Date(d)
	q=Dates.quarter(d)
	if freq=="Month"
		return(string(y)*"m"*string(m))
	elseif freq=="Quarter"
		return(string(y)*"Q"*string(m))
	elseif freq=="Year"
		return(string(y))
	end

end

MSA_CU=CSV.read("MSA.csv", DataFrame)
rent, cpi = download_price_index(MSA_CU)

d_rent=diff(rent, 1, true)
d_cpi=diff(cpi, 1, true)

m_rate = download_mortgage_rate()
d_m_rate=diff(m_rate, 1, true)


SP_housing = download_SP_housing_index()
SP_housing = @chain SP_housing begin
	@rsubset :Index_type == "XRSA"   #"XRNSA"
	@select(:date, :SP_area, :value)
end

wide_SP_housing=unstack(SP_housing, :date, :SP_area, :value)

dropmissing!(wide_SP_housing, disallowmissing = true)
d_housing=diff(wide_SP_housing, 1, true)

MSAs = unique(SP_housing[:, :SP_area])
df_MSA = DataFrame(SP_area = MSAs)
unemp_rate, unemp, emp, labor_force, area_unemp = download_unemployment(df_MSA)

d_unemp_rate=diff(unemp_rate, 1, true)
d_emp=diff(emp, 1, true)


# Emp = download_employment(area_unemp)

# long_temp = stack(temp, 2:ncol(temp))
# DataFramesMeta.rename!(long_temp, Symbol.(["date", "MSA", "index"]))

# long_var = @chain long_temp begin
# 	@transform(:Quarter = Dates.quarterofyear.(:date))
# 	@transform(:Year = Dates.year.(:date))
# 	@by([:MSA, :Year, :Quarter], :VAR = Statistics.var(:index), :stdev = Statistics.std(:index), :Q_m = mean(:index))
# 	@transform(:strSeq = string.(:Year) .* "Q" .* string.(:Quarter))
# 	@rtransform(:date = Dates.Date(:Year, (:Quarter - 1) * 3 + 1, 1))
# end

# wide_var = unstack(long_var, :date, :MSA, :Q_m)

# tables = DY12(wide_var)
# SpillPap.printTable(tables, names(wide_var)[2:end])

# CSV.write("./lib/spillover/var.csv", wide_var)

# num_MSAs = length(MSAs)

# average_length = Int(num_rows / num_MSAs) - 2
# wide_table_F = Matrix{Float64}(undef, average_length, num_MSAs)
# wide_table_B = Matrix{Float64}(undef, average_length, num_MSAs)
# log_price_to_rent = Matrix{Float64}(undef, Int(num_rows / num_MSAs), num_MSAs)
# seq = Vector{Float64}(undef, average_length)

bubbles=DataFrame[]
fundamentals=DataFrame[]

for index = eachindex(MSAs)

	msa::String = MSAs[index]
	println(msa)
	#df=DataFrame(date = Date[], housing = Float64[], rent=Float64[], cpi=Float64[], )
	#mat=
	#partial = df1[(df1.brief_area.==msa), [:y, :seq, :cpi, :log_emp, :une_rate, :log_rent, :log_population, :m_rate, :log_Close]]  #:log_num_permits
	df0=d_housing[:, [:date, Symbol(msa)]]
	DataFramesMeta.rename!(df0, [:date, :housing])


	df1=d_rent[:, [:date, Symbol(msa)]]
	DataFramesMeta.rename!(df1, [:date, :rent])

	df2=d_cpi[:, [:date, Symbol(msa)]]
	DataFramesMeta.rename!(df2, [:date, :cpi])
	
	df3=d_unemp_rate[:, [:date, Symbol(msa)]]
	DataFramesMeta.rename!(df3, [:date, :unemp])

	df4=d_emp[:, [:date, Symbol(msa)]]
	DataFramesMeta.rename!(df4, [:date, :emp])
	
	df=innerjoin(df0, df1, on=:date)
	df=innerjoin(df, df2, on=:date)
	df=innerjoin(df, df3, on=:date)
	df=innerjoin(df, df4, on=:date)	
	df=innerjoin(df, d_m_rate, on=:date)
	
	# num_obs, num_vars = size(partial)
	# if index == 1
	#   seq = Vector(partial.seq)[2:num_obs]
  
	# end
	temp_y = Vector(df.housing)  # 1 to num_obs  
	temp_x = Matrix(df[:, 2:end])
	dates=Vector(df.date)
	# log_price_to_rent[:, index] = temp_y
	AR = true
	#wide_table_F[:, index], wide_table_B[:, index] = _ivx.decompose(temp_y, temp_x, AR, 4, 1, 0.3, 0.02)
	println("-------")
	println(size(temp_x))
	F, B = _ivx.decompose(temp_y, temp_x, AR, 4, 1, 0.3, 0.02)

	F=DataFrame(date=dates[3:end], value=F)
	DataFramesMeta.rename!(F, [:date, Symbol(msa)])

	B=DataFrame(date=dates[3:end], value=B)
	DataFramesMeta.rename!(B, [:date, Symbol(msa)])

	push!(fundamentals, F)
	push!(bubbles, B)
	
  end

  df=bubbles[1]
  ## 时间有问题！！！！！！！！！！！！！！！！！！！！！！！！！！！
  for i in 2:length(bubbles)
	df=innerjoin(df, bubbles[i], on=:date)
  end




  y = df
  CSV.write( "newBubbles.csv",y)
  #wide_mean=unstack(long_var, strSeq, :MSA, :Q_m)
  strSeq = string.(y[:, 1])
  T, col_width = size(y)
  
  trend = false
  IC = "bic"
  
  max_lag = 4
  window_size = floor(Int64, 0.2 * T)  #18
  # 12              #% 12 years
  robust = true
  alpha = 0.95
  control_size = 6  #12
  
  results = _tvgc.tvgc(y, trend, max_lag, window_size, control_size, alpha, robust, IC)
  
  

#m=ARCHModels.fit(DCC{1, 1, GARCH{1, 1}}, data, method=:twostep)#, meanspec = NoIntercept)

#y=unstack(long_var, [:Year, :Quarter], :MSA, :Q_m)

num_plots = size(results[5], 1)

the_folder = "./png/TVGC/recursive/"
rm(the_folder, force = true, recursive = true)
mkdir(the_folder)

for i in 1:num_plots
	the_label = results[5][i, 1] * "\$\\longmapsto\$" * results[5][i, 2]
	plot(strSeq[window_size:T], results[1][:, i], label = the_label, dpi = 600, lw = 2, linecolor = :black)
	plot!(strSeq[window_size:T], results[2][:, i], label = "95% cv", line = (:dot, 2))
	savefig(the_folder * results[5][i, 1] * " to " * results[5][i, 2] * ".png")
end

the_folder = "./png/TVGC/rolling window/"
rm(the_folder, force = true, recursive = true)
mkdir(the_folder)

for i in 1:num_plots
	the_label = results[5][i, 1] * "\$\\longmapsto\$" * results[5][i, 2]
	plot(strSeq[window_size:T], results[3][:, i], label = the_label, dpi = 600, lw = 2, linecolor = :black)
	plot!(strSeq[window_size:T], results[4][:, i], label = "95% cv", line = (:dot, 2))
	savefig(the_folder * results[5][i, 1] * " to " * results[5][i, 2] * ".png")
end


bsadf, cv, d_label_original=gsadf.exuber(y,1,0)

col_width=size(y,2)
for i in 2:col_width

	d_label=format_date.(d_label_original, "Month")
    plot(d_label, bsadf[:,i-1], label="bsadf", dpi=1200, lw=2, linecolor=:black)
    plot!(d_label, cv[:,i-1], label="95% cv", line=(:dot, 2), linecolor=:red, dpi=1200)
    savefig("./png/Bubbles/" * names(y)[i] * ".png")

end




# #housing_index = CSV.read(download("https://www.freddiemac.com/fmac-resources/research/docs/fmhpi_master_file.csv"), DataFrame)



# con = DBInterface.connect(DuckDB.DB)
# # # create a DataFrame
# # # register it as a view in the database
# DuckDB.register_data_frame(con, housing_index, "my_df")
# CBSA_index = DBInterface.execute(con, "SELECT GEO_Code, Index_NSA, Index_SA, (Year-1975)*12+Month AS seq FROM my_df where GEO_Type ='CBSA' ")


# MSA_index = @chain housing_index begin
# 	@rsubset(:GEO_Type == "CBSA")
# 	@transform(:seq = (:Year .- 1975) * 12 + :Month)
# 	@transform(:seq_lag = :seq .+ 12)
# 	@select(:GEO_Code, :Year, :Month, :seq, :seq_lag, :Index_NSA, :Index_SA)
# end

# temp = innerjoin(MSA_index, MSA_index, on = [:seq => :seq_lag, :GEO_Code => :GEO_Code], renamecols = "_left" => "_right")


# growth = @chain temp begin
# 	@transform(:growth = log.(:Index_NSA_left) - log.(:Index_NSA_right))
# 	@transform(:Quarter = trunc.(Int32, (:Month_left .+ 2) / 3))
# 	@by([:GEO_Code, :Year_left, :Quarter], :VAR = Statistics.var(:growth))
# 	@select(:GEO_Code, :Year = :Year_left, :Quarter, :VAR)
# end

# @transform!(growth, :Quarter = trunc.(Int32, (:Month .+ 2) / 3))


# mortgage = fred("MORTGAGE30US")

# price_index = CSV.read(download("https://download.bls.gov/pub/time.series/cu/cu.data.0.Current"), DataFrame; normalizenames = true)
# @transform!(price_index, :area_code = SubString.(:series_id, 5, 8))
# @rtransform!(price_index, :item_code = SubString(:series_id, 9, length(:series_id)))
# @rtransform!(price_index, :seasonal = SubString(:series_id, 3, 3))
# @rtransform!(price_index, :period = SubString(:series_id, 4, 4))

# con = DBInterface.connect(DuckDB.DB)

# # create a DataFrame


# # register it as a view in the database
# DuckDB.register_data_frame(con, price_index, "my_df")

# # run a SQL query over the DataFrame
# results = DBInterface.execute(con, "SELECT * FROM my_df where area_code GLOB '0[1-9][1-9]0' and item_code LIKE 'SEHA%' ")

# subregions_rent = @chain price_index begin
# 	@rsubset contains(:area_code, r"0[1-9][1-9]0") contains(:item_code, r"SEHA")

# end


# subregions_rent = @chain price_index begin
# 	@rsubset contains(:area_code, r"0[1-9]00") contains(:item_code, r"SEHA")
# end

# s_regions_rent = @chain price_index begin
# 	@rsubset contains(:area_code, r"[S][1-9]00") contains(:item_code, r"SEHA")
# end
