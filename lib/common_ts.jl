module _common_TS

using LinearAlgebra
using Statistics
using Random
using HypothesisTests

export regress, ADF, get_reg

function regress(x, y)
#y only consists of one column
    dof = size(y, 1) - size(x, 2)
    x_t = x'
    temp = (x_t * x)
    g = inv(temp)
    beta = g * (x_t * y)                               # @-model A-@
    eps = y - x * beta

    se = eps' * eps / dof

    sig = sqrt.(diag(se * g))

    return (beta, eps, se, sig)
end

function get_reg(y, deterministic::Symbol, lag)#, lag_first=true)
    obs, num_y = size(y)  #bug: must be a matrix
    #return的矩阵，实际从lag+1开始
    if deterministic==:none
        num_x = 0
        reg = Matrix{Float64}(undef, obs, num_y * lag)    
    elseif deterministic == :constant
        num_x = 1
        reg = [ones(obs) Array{Float64}(undef, obs, num_y * lag)]
    elseif deterministic == :trend
        num_x = 2
        reg = [ones(obs) @view(x[1:obs]) Array{Float64}(undef, obs, num_y * lag)]
    else
        throw(ArgumentError("deterministic = $(deterministic) is invalid"))
    end
   
    @inbounds @simd for j ∈ 1:lag
        reg[j+1:obs, num_x+num_y*(j-1)+1:num_x+num_y*j] = @view(y[1:obs-j, :])
    end
    
    return reg
end

function ADF(y, IC, deterministic::Symbol, adflag, null_hyp=false)

    T0 = size(y,1)
        
    dy = diff(y) #@view(y[2:T0]) - @view(y[1:T0-1]) #y(1:T1)    

    real_T = T0 - 1- adflag   #from adflag+1 to T0-1

    reg_full = get_reg(@view(dy[:, 1:1]), deterministic , adflag)

    lag_y =@view y[adflag+1:T0-1, :]         #@-from k+1 to the end (including y1 and x)-@
    dy01 =@view  dy[adflag+1:end]      #@-from k+1 to the end (including dy0)-@

    if IC =="bic" || IC=="aic"   #  IC: 0 for fixed lag order 1 for AIC and 2 for BIC
        ICC = zeros(adflag + 1)
        ADF = zeros(adflag + 1)
        beta_list = Any[]  #cell(adflag+1,1);
        eps_list = Any[]  #cell(adflag+1,1);

        for k = 0:adflag
            # Model specification 
            if null_hyp
                x = @view(reg_full[adflag+1:end, 1:(k+1)])
            else
                x = [lag_y @view(reg_full[adflag+1:end, 1:(k+1)])]
            end

            # OLS regression
            beta, eps, se, sig = regress(x, dy01)
            push!(beta_list, beta)
            push!(eps_list, eps)

            # Information Criteria
            npdf = sum(-1 / 2 * log(2 * pi) .- 1 / 2 * (eps .^ 2))
            if IC == "aic"               #@ AIC @
                ICC[k+1] = -2 * npdf / real_T + 2 * size(beta, 1) / real_T
            elseif IC == "bic"           #@ BIC @
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
            x = @view(reg_full[adflag+1:end, 1:(adflag+1)])
        else
            x = [lag_y @view(reg_full[adflag+1:end, 1:(adflag+1)])]
        end
        beta, eps, se, sig = regress(x, dy01)
        ## OLS regression
        ADFlag = beta[1] / sig[1]

        return (ADFlag, beta, eps, adflag)

    end


end
end