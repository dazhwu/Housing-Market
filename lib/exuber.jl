
__precompile__(true)
module gsadf
using DataFrames
using LinearAlgebra
using Statistics
using Random
using DelimitedFiles
using Dates


export rls_gsadf, CV_PSY, gsadf_results, wmboot, exuber

struct gsadf_results
	num_points::Int64
	start_point::Int64
	end_point::Int64  #ending period
	badf::Vector{Float64}
	bsadf::Vector{Float64}
	adf::Float64
	sadf::Float64
	gsadf::Float64
end

function exuber(df::DataFrame, adflag::Int64, Min_Window::Int64=0)
	obs = nrow(df)
	if Min_Window<=0
		r0 = 0.01 + 1.8 / sqrt(obs)
		Min_Window = floor(Int64,r0 * obs)
	end
	
	d_label = df[:,1][Min_Window+adflag+1:obs]
	bsadf=Matrix{Float64}(undef, obs-Min_Window-adflag, ncol(df)-1)
	cv=Matrix{Float64}(undef, obs-Min_Window-adflag, ncol(df)-1)
#	d_label_old = strSeq[swindow0:T]
	for i in 2:ncol(df)
		ts=df[:,i]
		results = rls_gsadf(ts, adflag, Min_Window)
		bsadf[:,i-1]=results.bsadf
    	cv[:,i-1] =@view CV_PSY(obs, 1000, Min_Window, adflag)[:,2]

	end
	return (bsadf, cv, d_label)
end

function get_reg!(y::Vector{T}, lag_y::Vector{T}, lag::Int64, reg::Matrix{T}, null_hyp=false) where {T <: Float64}
#each row is a variable and each col is an observation
	ncol = size(reg, 2)
	start_row=ifelse(null_hyp, 1,2)

	@inbounds @simd for i in 1:ncol  #eachcol(reg)    # eachc1:ncol   #column
		
		if null_hyp==false
			reg[1, i] = lag_y[i]	
		end
		reg[start_row, i] = 1
			@inbounds @simd for j in 1:lag
				reg[j+start_row, i] = y[lag+1+i-j] - y[lag+i-j]
			end
		end
end

function rls_zero(y::Vector{T}, x::Vector{T}, tstat::Matrix{T}, max_start::Int64, min_win::Int64, obs::Int64) where {T <: Float64}


	x_x = Vector{T}(undef, length(x))
	x_y = Vector{T}(undef, length(x))
	@inbounds @simd for i in eachindex(x)
		x_x[i] = x[i] * x[i]
		x_y[i] = x[i] * y[i]
	end
	@inbounds @simd for start ∈ 1:max_start
		sumx = zero(T)
		sumy = zero(T)
		sumxx = zero(T)
		sumxy = zero(T)
		@inbounds for win_size ∈ min_win:obs-start+1
			end_point = start + win_size - 1
			if win_size == min_win #&& start==1

				@inbounds for i in start:end_point
					sumx += x[i]
					sumy += y[i]
					sumxx += x_x[i]
					sumxy += x_y[i]
				end

			else
				sumx += x[end_point]
				sumy += y[end_point]
				sumxx += x_x[end_point] #* x[end_point]
				sumxy += x_y[end_point]  #* x[end_point]

			end
			meanx = sumx / win_size
			meany = sumy / win_size
			den = sumxx / win_size - meanx * meanx
			beta = (sumxy / win_size - meanx * meany) / den
			alpha = meany - beta * meanx

			u = Vector{T}(undef, end_point - start + 1)
			@inbounds @simd for i in eachindex(u)
				u[i] = y[start+i-1] - alpha - beta * x[start+i-1]
			end

			suu = u'u
			sbeta = sqrt(suu / (win_size - 2) / den / win_size)
			tstat[start, end_point-min_win+1] = beta / sbeta
		end
	end

end
function rls_multi(y, t_x::Matrix{T}, tstat::Matrix{T}, max_start::Int64, min_win::Int64, obs::Int64) where {T <: Float64}

	x_width = size(t_x, 1)

	#t_x = x'

	g = Matrix{T}(undef, x_width, x_width)
	b = Vector{T}(undef, x_width)


	@inbounds @simd for start ∈ 1:max_start
		@inbounds for win_size ∈ min_win:obs-start+1
			end_point = start + win_size - 1

			tsx = @view t_x[:, start:end_point]

			sy = @view y[start:end_point]

			if win_size == min_win

				g = inv(tsx * tsx')
				b = g * (tsx * sy)
			else
				t_new_x = t_x[:, [end_point]]
				new_x = t_new_x'

				gx = g * t_new_x
				the_factor = 1 / (1 + (new_x*gx)[1])
				g -= the_factor * (gx * gx')

				temp = -y[end_point]
				@inbounds @simd for i in eachindex(b)
					temp += b[i] * t_new_x[i]
				end
				#b -= g*t_new_x*(dot(b, t_new_x)-y[end_point])
				b -= g * t_new_x * temp
			end
			sqres = zero(T)
			res = Vector{T}(undef, win_size)
			@inbounds @simd for i in eachindex(sy)
				temp = zero(T)
				@inbounds @simd for j in eachindex(b)
					temp += b[j] * tsx[j, i]
				end
				res[i] = sy[i] - temp #dot(tsx[:,i], b)
				sqres += res[i] * res[i]

			end
			# res = sy - tsx'b
			# sqres = (res'res)[1]

			vares::T = sqres / (win_size - x_width)
			sb_1::T = sqrt(vares * diag(g)[1])


			tstat[start, end_point-min_win+1] = (b[1]) / sb_1
			#    write(file, tstat[i,j])
		end
	end

end


function rls_gsadf(y::Vector{T}, adflag::Int64, min_win::Int64) where {T <: Float64}

	obs_ori = length(y)
	if min_win == 0
		r0 = 0.01 + 1.8 / sqrt(obs_ori)
		min_win = floor(Int64, r0 * obs_ori)
	end

	dy01 = Vector{T}(undef, obs_ori - 1 - adflag)
	lag_y = Vector{T}(undef, obs_ori - 1 - adflag)
	@inbounds @simd for i in 1:obs_ori-1-adflag

		lag_y[i] = y[adflag+i]

		dy01[i] = y[adflag+1+i] - y[adflag+i]

	end

	#dy01 = @view dy[adflag+1:end]
	obs = length(dy01)

	max_start::Int64 = obs - min_win + 1  #max_starting

	tstat = Matrix{T}(undef, max_start, max_start)

	if adflag == 0
		rls_zero(dy01, lag_y, tstat, max_start, min_win, obs)

	else
		x = Matrix{T}(undef, 2 + adflag, obs_ori - adflag - 1)

		get_reg!(y, lag_y, adflag, x)
		rls_multi(dy01, x, tstat, max_start, min_win, obs)
	end




	adf = tstat[1, max_start][1]
	badf = @view(tstat[1, :])
	sadf = findmax(badf)[1]

	bsadf = Vector{T}(undef, max_start)
	@inbounds @simd for i in 1:max_start
		bsadf[i] = findmax(@view(tstat[1:i, i]))[1]
	end


	gsadf = findmax(bsadf)[1]

	return (gsadf_results(max_start, adflag + 2, obs_ori, badf, bsadf, adf, sadf, gsadf))



end


function CV_PSY(T::Int, m::Int, swindow0::Int, adflag::Int)

	qe = [0.90 0.95 0.99]


	dim = T - adflag - swindow0
	MPSY = Matrix{Float64}(undef, dim, m)
	tbr = Matrix{Float64}(undef, dim, length(qe))


	seed = 123


	Random.seed!(seed)
	e = randn(T, m)
	a = T^(-1)

	y = cumsum(e .+ a, dims = 1)


	@inbounds @simd for iter ∈ 1:m
		if mod(iter, 100) == 0
			println(iter)
		end
		result = rls_gsadf(y[:, iter], adflag, swindow0)             #PSY_fixed(y[:, iter], swindow0, adflag)
		MPSY[:, iter] = accumulate(max, result.badf)
	end


	@inbounds @simd for i in 1:dim
		@inbounds for j in eachindex(qe)

			tbr[i, j] = quantile(@view(MPSY[i, :]), qe[j])
		end
	end
	return (tbr)
end

function regress(x, y)

	dof = size(y, 1) - size(x,1)
	
	temp = (x * x')
       g=inv(temp)
	beta = g * (x * y)                               # @-model A-@
	eps = y - x' * beta

	#se = eps' * eps / dof

	#sig = sqrt.(diag(se * g))

	return (beta, eps)
end


function wmboot(y, min_win, adflag, nboot=499)

	qe = [0.90 0.95 0.99]
	
	obs_ori = length(y)

	dim = obs_ori - adflag - min_win
	MPSY = Matrix{Float64}(undef, dim, nboot)
	tbr = Matrix{Float64}(undef, dim, length(qe))

	dy01 = Vector{Float64}(undef, obs_ori - 1 - adflag)
	
	@inbounds @simd for i in 1:obs_ori-1-adflag

		

		dy01[i] = y[adflag+1+i] - y[adflag+i]

	end

	#dy01 = @view dy[adflag+1:end]
	obs = length(dy01)
	
	reg = Matrix{Float64}(undef, 1 + adflag, obs_ori - adflag - 1)
	get_reg!(y, Vector{Float64}(undef,1), adflag, reg, true)
	



	beta, eps = regress(reg, dy01)
	lag=adflag
	T0 = length(eps)

	T = size(y, 1)
	dy = y[2:T] - y[1:T-1]

	Tb=T
	Random.seed!(6)
	rN = rand(1:T0, Tb, nboot)
	wn = randn(Tb, nboot)

	dyb = zeros(Tb - 1, nboot)
	dyb[1:lag, :] = repeat(dy[1:lag], 1, nboot)

	for j = 1:1:nboot
		   if lag == 0
				  for i = 1:Tb-1
						 dyb[i, j] = wn[i, j] * eps[rN[i, j]]
				  end

		   elseif lag > 0
				  x = zeros(Tb - 1, lag)
				  for i = lag+1:1:Tb-1
						 for k = 1:1:lag
								x[i, k] = dyb[i-k, j]
						 end
						
						 dyb[i, j] = dot(x[i, :], beta[2:end]) + wn[i-lag, j] * eps[rN[i-lag, j]]
				  end
		   end
	end

	#yb0 = ones(1, nboot)*y[1]
	dyb0 = vcat(ones(1, nboot) * y[1], dyb)
	yb = cumsum(dyb0, dims=1)

	@inbounds @simd for iter ∈ 1:nboot
		if mod(iter, 100) == 0
			println(iter)
		end
		result = rls_gsadf(yb[:, iter], adflag, min_win)             #PSY_fixed(y[:, iter], swindow0, adflag)
		MPSY[:, iter] = accumulate(max, result.badf)
	end


	@inbounds @simd for i in 1:dim
		@inbounds for j in eachindex(qe)

			tbr[i, j] = quantile(@view(MPSY[i, :]), qe[j])
		end
	end
	return (tbr)

	## ####The PSY Test########################


	# ##parpool(6);
	# dim = Tb - swindow0 + 1
	# MPSY = zeros(dim, nboot)
	# for iter = 1:nboot
	# 	   MPSY[:, iter] = PSY(@view(yb[:, iter]), swindow0, IC, adflag)
	# end
	# println("wmboot_psy ", Dates.now() - start)
	# ####delete(gcp);

	# SPSY = vec(findmax(MPSY, dims=1)[1])  #max of each col
	# tbr=quantile(SPSY, qe)
	# println(Dates.now() - start)
	# return (tbr)

end





end
