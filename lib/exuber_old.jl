module exuber_old
using LinearAlgebra
using Statistics
using Random
using DelimitedFiles
using Dates
#include("tvgc.jl")


export ADF,wmboot



function regress(x, y)

	dof = size(y, 1) - size(x,2)
	x_t = x'
	temp = (x_t * x)
       g=inv(temp)
	beta = g * (x_t * y)                               # @-model A-@
	eps = y - x * beta

	se = eps' * eps / dof

	sig = sqrt.(diag(se * g))

	return (beta, eps, se, sig)
end


function ADF(y, IC, adflag, null_hyp=false)

       T0 = length(y)
       T1 = T0 - 1    #T1=size(y,1)-1   #T0 T1 size(y,1)

       y1 = @view(y[1:T1])
       dy = @view(y[2:T0]) - y1 #y(1:T1)
       x1 = [y1 ones(T1, 1)]

       real_T = T1 - adflag

       reg_full = _tvgc.get_reg(@view(dy[:, 1:1]), [0], adflag)

       xx = y[adflag+1:T1, :]         #@-from k+1 to the end (including y1 and x)-@
       dy01 = dy[adflag+1:T1]      #@-from k+1 to the end (including dy0)-@

       if IC > 0   #  IC: 0 for fixed lag order 1 for AIC and 2 for BIC
              ICC = zeros(adflag + 1)
              ADF = zeros(adflag + 1)
              beta_list = Any[]  #cell(adflag+1,1);
              eps_list = Any[]  #cell(adflag+1,1);

              for k = 0:1:adflag
                     # Model specification 
                     if null_hyp
                            x2 = @view(reg_full[adflag+1:end, 1:(k+1)])
                     else
                            x2 = [xx @view(reg_full[adflag+1:end, 1:(k+1)])]
                     end

                     # OLS regression
                     beta, eps, se, sig = regress(x2, dy01)
                     push!(beta_list, beta)
                     push!(eps_list, eps)

                     # Information Criteria
                     npdf = sum(-1 / 2 * log(2 * pi) .- 1 / 2 * (eps .^ 2))
                     if IC == 1               #@ AIC @
                            ICC[k+1] = -2 * npdf / real_T + 2 * size(beta, 1) / real_T
                     elseif IC == 2           #@ BIC @
                            ICC[k+1] = -2 * npdf / real_T + size(beta, 1) * log(real_T) / real_T
                     end



                     ADF[k+1] = beta[1] / sig[1]
              end

              _, lag0 = findmin(ICC)

              ADFlag = ADF[lag0]

              return (ADF[lag0], beta_list[lag0], eps_list[lag0], (lag0 - 1))   #because it starts with 0
              #lag=lag0-1

       else
              if null_hyp
                     x2 = @view(reg_full[adflag+1:end, 1:(adflag+1)])
              else
                     x2 = [xx @view(reg_full[adflag+1:end, 1:(adflag+1)])]
              end
              beta, eps, se, sig = regress(x2, dy01)
              ## OLS regression


              ADFlag = beta[1] / sig[1]

              return (ADFlag, beta, eps, adflag)

       end

end

function PSY(y, swindow0, IC, adflag)

       T = length(y)
       bsadfs = Matrix{Float64}(undef, T - swindow0 + 1, 1)
       for r2 = swindow0:T
              rwadft = zeros(r2 - swindow0 + 1, 1)
              for r1 = 1:r2-swindow0+1
                     rwadft[r1], _, _, _ = ADF(y[r1:r2, 1:1], IC, adflag)   ## two tail 5# significant level
              end
              bsadfs[r2-swindow0+1] = findmax(rwadft)[1]
       end

       return (bsadfs)
end


function CV_PSY(T, swindow0, IC, adflag)

       qe = [0.90; 0.95; 0.99]
       m = 100

       dim = T - swindow0 + 1
       MPSY = zeros(dim,m)

       SI = 1
       start = Dates.now()

       Random.seed!(SI)
       e = randn(T, m)
       a = T^(-1)
       y = cumsum(e .+ a, dims=1)
       println(Dates.now() - start)
       @inbounds @simd  for iter = 1:m
              println(iter)
              MPSY[:, iter] = PSY(y[:, iter], swindow0, IC, adflag)
       end
       println("cv_psy", Dates.now() - start)
       # SPSY = vec(findmax(MPSY, dims=1)[1])  #max of each col
       # tbr=quantile(SPSY, qe)
       return ([quantile(@view(MPSY[i,:]), qe) for i in 1:size(MPSY,1)])
end



function wmboot(y, swindow0, IC, adflag, Tb, nboot=499)

       qe = [0.90; 0.95; 0.99]

       _, beta, eps, lag = ADF(y, IC, adflag, true)
       T0 = length(eps)

       T = size(y, 1)
       dy = y[2:T] - y[1:T-1]

       ## The DGP
       start = Dates.now()
       Random.seed!(6)
       rN = rand(1:T0, Tb, nboot)
       wn = randn(Tb, nboot)

       dyb = zeros(Tb - 1, nboot)
       dyb[1:lag, :] = repeat(dy[1:lag], 1, nboot)

       for j = 1:1:nboot
              if lag == 0
                     for i = lag+1:1:Tb-1
                            dyb[i, j] = wn[i-lag, j] * eps[rN[i-lag, j]]
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


       println(Dates.now() - start)
       ## ####The PSY Test########################


       ##parpool(6);
       dim = Tb - swindow0 + 1
       MPSY = zeros(dim, nboot)
       for iter = 1:nboot
              MPSY[:, iter] = PSY(@view(yb[:, iter]), swindow0, IC, adflag)
       end
       println("wmboot_psy ", Dates.now() - start)
       ####delete(gcp);

       SPSY = vec(findmax(MPSY, dims=1)[1])  #max of each col
       tbr=quantile(SPSY, qe)
       println(Dates.now() - start)
       return (tbr)

end



end