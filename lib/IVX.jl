

__precompile__(true)

module _ivx

#using ShiftedArrays
using Statistics
using LinearAlgebra
using Distributions
using StateSpaceModels
using HypothesisTests

export IVX, IVX_AR, ivx_ar, ivx, decompose

mutable struct IVX
  aivx::AbstractArray{Float64}
  wivxind::Matrix{Float64}
  p_value_ivx::Matrix{Float64}
  wivxind_z::Matrix{Float64}
  df::Int64
  horizons::Int64
  df_residuals::Int64
  delta::Matrix{Float64}
  varcov::Matrix{Float64}
  Rn::Vector{Float64}
  Rz::Vector{Float64}
  intercept::Matrix{Float64}
  fitted::Matrix{Float64}
  residuals::Matrix{Float64}
  interceptm::Matrix{Float64}
  fittedm::Matrix{Float64}
  ols_residuals::Vector{Float64}

end



mutable struct IVX_AR
  ivx::IVX   #Ref{IVX}
  order_p::Int64
  rse::Float64
  coefficients_ar::Vector{Float64}
end


function ivx_ar(y::AbstractArray{Float64}, x::AbstractArray{Float64}, ar_max, Horizon=1, offset=0.3, step=0.02)

  mdl_ivx = ivx(y, x, Horizon)
  mdl_ar = StateSpaceModels.auto_arima(mdl_ivx.ols_residuals, max_p=5, max_q=0, max_d=0, d=0, show_trace=false, information_criteria="bic")

  order_p = mdl_ar.order.p

  ar_coefs = mdl_ar.results.coef_table.coef[1:order_p]
  if (order_p == 0)
    return (IVX_AR(mdl_ivx, 0, 0, [0])) #(res_ivx[chosen])
  end
  ngrid::Int64 = offset / step * 2 + 1

  grid_seq = Array{Float64}(undef, ngrid, order_p)

  for j = 1:order_p
    @inbounds @simd for i = 1:ngrid
      grid_seq[i, j] = ar_coefs[j] - offset + (i - 1) * step
    end

  end
  #println(grid_seq)
  #ngrid, _ = size(grid_seq)
  min_rse = typemax(Float64)
  res_ivx = Vector{IVX}(undef, ngrid)
  #res_ivx=Ref{IVX}
  chosen::Int64 = -1
  for i = 1:ngrid

    x_adj = tilt(x, grid_seq[i, :], order_p)
    y_adj = tilt(@view(y[:, 1:1]), grid_seq[i, :], order_p)

    res_ivx[i] = ivx(y_adj, x_adj, Horizon)

    eps = y_adj .- sum(x_adj * res_ivx[i].aivx)
    rse = Statistics.var(eps)

    if rse < min_rse
      min_rse = rse
      chosen = i
    end
  end

  return (IVX_AR(res_ivx[chosen], order_p, min_rse, grid_seq[chosen, :])) #(res_ivx[chosen])

  # z$Wald_AR <- ac_test_wald(mdl_ivx$ols$residuals, q)
  # z$q <- q
  # z


end

function ivx(y::AbstractArray{Float64}, x::AbstractArray{Float64}, K=1)

  nr = length(y) #size(y,1) #number of rows

  xlag = x[1:end-1, :]  #height: nr-1

  xt = x[2:end, :]
  yt = y[2:end, :]
  num_x_var = size(x, 2)       #number of x
  nn = nr - 1
  Xols = hcat(ones(nn), xlag)

  Aols = (Xols' * Xols) \ (Xols' * yt)                #(Xols \ yt)
  epshat = yt - Xols * Aols
  s2 = sum(epshat .^ 2, dims=1)

  std_err = sqrt.(s2 .* diag(inv(Xols' * Xols)))
  tstat = Aols ./ std_err

  Rn = zeros(Float64, num_x_var, num_x_var)
  for i in 1:num_x_var
    coef_Rn = (xlag[:, i] \ xt[:, i])    #vector on vector #[1,1]           #as_scalar(arma::solve(xlag.col(i), xt.col(i)));        
    Rn[i, i] = coef_Rn
  end
  u = xt - xlag * Rn  #number of columns in u: L    
  corrmat = Statistics.cor(epshat, u)

  covepshat = epshat' * epshat / nn   #result: 1x1 matrix

  covu = zeros(Float64, num_x_var, num_x_var)

  covu = (u' * u) / nn

  covuhat = u' * epshat / nn

  m = floor(Int64, nn^(1 / 3))
  uu = zeros(num_x_var, num_x_var)
  for h = 1:m
    u_h = u[(h+1):nn, :]
    u_0 = u[1:(nn-h), :]
    a = u_0' * u_h
    uu += (1 - h / (m + 1)) * a
  end
  uu = uu / nn
  Omegauu = covu + uu + transpose(uu)
  q = zeros(1, num_x_var)
  for h in 1:m
    p = zeros(1, num_x_var)
    for t = 1:(nn-h)
      p[1, :] += u[t+h, :] * epshat[t]
    end

    q += ((1 - h / (m + 1)) * p)
  end
  residue = q / nn
  Omegaeu = covuhat + residue'
  Rz = (1 - 1 / (nn^0.95)) * Matrix{Float64}(I, num_x_var, num_x_var)
  diffx = xt - xlag
  z = zeros(nn, num_x_var)

  z[1:1, :] = diffx[1:1, :]

  for i in 2:nn

    z[[i], :] = z[[i - 1], :] * Rz + diffx[[i], :]
  end

  n = nn - K + 1
  Z = vcat(zeros(1, num_x_var), z[1:(n-1), :])
  zz = vcat(zeros(1, num_x_var), z[1:(nn-1), :])
  ZK = zeros(n, num_x_var)
  for i = 1:n
    ZK[i:i, :] = sum(zz[i:(i+K-1), :], dims=1)
  end

  meanzK = Statistics.mean(ZK, dims=1)
  #######################
  yy = zeros(n, 1)
  for i = 1:n
    yy[i:i, :] = sum(@view(yt[i:(i+K-1), :]), dims=1)
  end
  Yt = yy .- mean(yy, dims=1)
  xK = zeros(n, num_x_var)
  for i = 1:n
    xK[i:i, :] = sum(@view(xlag[i:(i+K-1), :]), dims=1)
  end
  meanxK = Statistics.mean(xK, dims=1)
  Xt = zeros(n, num_x_var)
  #print(size(ones(n,1)*meanxK[:,1]))

  for i in 1:num_x_var
    Xt[:, i] = xK[:, i] - ones(n, 1) * meanxK[:, i]
  end

  Aivx = Yt' * Z * inv(Xt' * Z)
  Aivx_t = Aivx'
  fitted = Xt * Aivx_t


  intercept = mean(Yt, dims=1) - mean(Xt, dims=1) * Aivx_t
  residuals = Yt - fitted
  ###########///////////////// No demeaning /////////////////
  interceptm = mean(y, dims=1) - mean(xlag, dims=1) * Aivx_t
  fittedm = interceptm .+ xlag * Aivx_t
  FM = covepshat - Omegaeu' * inv(Omegauu) * Omegaeu

  M = ZK' * ZK * (covepshat[1, 1]) - n * meanzK' * meanzK * FM[1, 1]
  H = Matrix{Float64}(I, num_x_var, num_x_var)
  Q = H * inv(Z' * Xt) * M * inv(Xt' * Z) * H'
  wivx = (H * Aivx_t)' * inv(Q) * (H * Aivx_t)

  wivxind_z = Aivx_t ./ (sqrt.(diag(Q)))
  wivxind = wivxind_z .^ 2
  p_value_ivx = ccdf(Chisq(1), wivxind)


  tbr = IVX(Aivx_t, wivxind, p_value_ivx, wivxind_z, num_x_var, K, nn - 1, corrmat, Q, diag(Rn), diag(Rz), intercept, fitted, residuals, interceptm, fittedm, @view(epshat[:, 1]))


  return tbr
end

function tilt(x::AbstractArray{Float64}, grid_vec::Array{Float64}, ar_length::Int64)

  nrow, ncol = size(x)

  tbr = Array{Float64}(undef, nrow - ar_length, ncol)
  temp_x = Array{Float64}(undef, nrow - ar_length, ar_length)
  for i in 1:ncol
    temp_y = @view(x[ar_length+1:nrow, i])
    for j in 1:ar_length
      temp_x[:, j] = @view(x[ar_length+1-j:nrow-j, i])
    end
    tbr[:, i] = temp_y - temp_x * grid_vec
  end

  return (tbr)
end


function decompose(y::AbstractArray{Float64}, predictor::AbstractArray{Float64},   AR::Bool, ar_max, horizon=1, offset=0.3, step=0.02)
  
  num_obs, num_vars=size(predictor)  

  growth=@view(y[2:num_obs])-@view(y[1:num_obs-1])
  
  x=@view(predictor[2:num_obs, :]) 

  coef=Vector{Float64}(undef, num_vars)
  if AR
    temp_ivx=_ivx.ivx_ar(growth, x, ar_max, horizon, offset, step)
    #ar_max, Horizon=1, offset=0.3, step=0.02
    coef=temp_ivx.ivx.aivx
    
  else
    temp_ivx=_ivx.ivx(growth, x)
    coef=temp_ivx.aivx
  end
  
  #temp=x*coef
  
  real_x=@view(predictor[2:num_obs-1, :])               #@view(x[1:num_obs-2,:])
  real_growth=@view(growth[2:end]) #num_obs-1])   start from y[3]-y[2]
  
  fitted_main=real_x * coef 
  intercept=mean(real_growth)-mean(fitted_main)
  
  fitted=fitted_main .+ intercept
  fundamental=Vector{Float64}(undef, num_obs-2)
  bubble=Vector{Float64}(undef, num_obs-2)
  @inbounds for i in eachindex(bubble)
    
    if i==1
      fundamental[i]=y[2]+fitted[1]
    else
      fundamental[i]=fundamental[i-1]+fitted[i]
    end
  end
  bubble=@view(y[3:num_obs])-fundamental
  
  return(fundamental, bubble)
end


end