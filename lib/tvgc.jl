module _tvgc

using DataFrames
using LinearAlgebra
using Statistics
using Random
using HypothesisTests
#using LoopVectorization
#include("exuber_old.jl")
include("common_ts.jl")

export tvgc, testSettings,  pair_tvgc

struct testSettings
    window_size::Int64
    control_size::Int64
    max_lag::Int64
    n_boot::Int64
    d::Int64
    alpha::Float32
    fast::Bool
    trend::Bool
    robust::Bool
end


function get_orders(y, trend)
    num_MSAs = size(y, 2)

    lag_orders = Vector{Int64}(undef, num_MSAs)
    println(size(lag_orders))

    for i in 1:num_MSAs
        println(i)
        lag_orders[i] = 0
        ts = y[:, i]
        rejected = false


        while rejected == false
            
            
            if trend
                _, _, _, the_lag =_common_TS.ADF(ts, "bic", :trend, 4)
                temp = (ADFTest(ts, :trend, the_lag))
            else
                _, _, _, the_lag =_common_TS.ADF(ts, "bic", :constant, 4)
                temp = (ADFTest(ts, :constant, the_lag))
            end
            
            if temp.stat < temp.cv[2]
                println("rejected")
                rejected = true
            else
                ts = diff(ts)  #@views ts[2:end] - ts[1:end-1]
                lag_orders[i] += 1
                println("+1")
            end
        end
    end


    return (lag_orders)

end


function tvgc(df, trend::Bool, max_lag::Int64, window_size, control_size, alpha, robust::Bool, IC)

    y = Matrix{Float64}(df[:, 2:ncol(df)])
    MSAs = names(df)[2:end]

    strSeq = Vector(df[:, 1])
    obs, num_MSAs = size(y)

    lag_order = get_orders(y, trend)
    result_rrGC = Matrix{Float64}(undef, obs - window_size + 1, num_MSAs * (num_MSAs - 1))

    rusult_rrCV = Matrix{Float64}(undef, obs - window_size + 1, num_MSAs * (num_MSAs - 1))

    result_rGC = Matrix{Float64}(undef, obs - window_size + 1, num_MSAs * (num_MSAs - 1))
    rusult_rCV = Matrix{Float64}(undef, obs - window_size + 1, num_MSAs * (num_MSAs - 1))
    labels = Matrix{String}(undef, num_MSAs * (num_MSAs - 1), 2)
    the_index = 1
    for i in 1:num_MSAs-1
        for j in i+1:num_MSAs

            d = max(lag_order[i], lag_order[j])
            temp_y = @view(y[:, [i, j]])
            result_1, result_2, header = pair_tvgc(temp_y, trend, max_lag, d, window_size, control_size, [[1, 2], [2, 1]], alpha, robust, IC)

            result_rrGC[:, the_index] = result_1[:, 9]
            rusult_rrCV[:, the_index] = result_1[:, 12]
            result_rrGC[:, the_index+1] = result_2[:, 9]
            rusult_rrCV[:, the_index+1] = result_2[:, 12]

            result_rGC[:, the_index] = result_1[:, 5]
            rusult_rCV[:, the_index] = result_1[:, 8]
            result_rGC[:, the_index+1] = result_2[:, 5]
            rusult_rCV[:, the_index+1] = result_2[:, 8]

            labels[the_index, 1] = MSAs[i]
            labels[the_index, 2] = MSAs[j]
            labels[the_index+1, 1] = MSAs[j]
            labels[the_index+1, 2] = MSAs[i]

            the_index += 2

        end
    end

    return (result_rrGC, rusult_rrCV, result_rGC, rusult_rCV, labels)

end


function pair_tvgc(y, trend::Bool, max_lag::Int64, d::Int64, window_size, control_size, FROM_TOs, alpha, robust::Bool, IC)
    settings::testSettings = testSettings(window_size, control_size, max_lag, 499, d, alpha, false, trend, robust)
    T = size(y, 1)



    n_windows = T - window_size + 1   #the number of windows
    M = window_size + control_size - 1
    if trend
        lag=var_soc(y, max_lag, :trend, IC)
        reg_full = _common_TS.get_reg(y, :trend, lag + d)
    else
        lag=var_soc(y, max_lag, :constant, IC)
        reg_full = _common_TS.get_reg(y, :constant, lag + d)
    end
    
    list_results = Any[]#Vector{Array{Float64}(dim, 12)}

    for i in eachindex(FROM_TOs)
        FROM = FROM_TOs[i][1]
        TO = FROM_TOs[i][2]
        push!(list_results, test_FROM_TO(y, reg_full, T, n_windows, M, lag, settings, FROM, TO))
    end




    # df = DataFrame(result_mat, :auto)
    # rename!(df, ["fGC", "fSC", "fPV", "fCV", "rGC", "rSC", "rPV", "rCV", "rrGC", "rrSC", "rrPV", "rrCV"])
    #rename!(df, ["fGC", "fSC", "fcv", "rGC", "rSC", "rcv", "rrGC", "rrSC", "rrcv"])
    #header="fGC, fSC, fPV , fCV, rGC, rSC, rPV, rCV, rrGC, rrSC, rrPV, rrCV"
    header = ["fGC", "fSC", "fPV", "fCV", "rGC", "rSC", "rPV", "rCV", "rrGC", "rrSC", "rrPV", "rrCV"]
    return list_results[1], list_results[2], header

end

function test_FROM_TO(y,  reg_full, nob, dim, M, lag, settings::testSettings, From_var, To_var)
    result_mat = Array{Float64}(undef, dim, 12)

    rwGC = Array{Float64}(undef, dim)
    rwSC = Array{Float64}(undef, dim)

    @inbounds for r2 in settings.window_size:nob
        dim0 = r2 - settings.window_size + 1   #dim0 will be 1, 2, 3, ..., dim

        @inbounds for r1 in 1:dim0

            temp_a, temp_b = granger_cause_Mwald(y, r1, r2, lag, From_var, To_var, reg_full, settings)
            rwGC[r1] = temp_a
            rwSC[r1] = temp_b
            if r1 == 1
                result_mat[dim0, 1] = temp_a
                result_mat[dim0, 2] = temp_b
            end
            if r1 == dim0   #cannot be combined with the previous if
                result_mat[dim0, 5] = temp_a
                result_mat[dim0, 6] = temp_b
            end
        end
        result_mat[r2-settings.window_size+1, 9], ind = findmax(@view(rwGC[1:dim0]))
        result_mat[r2-settings.window_size+1, 10] = rwSC[ind]
    end


    Mfcv, Mrcv, Mrrcv = bootstrap_RGC_MW(y, lag,  From_var, To_var, reg_full, M, settings)


    i = 1
    for max_element_each_row in [Mfcv, Mrcv, Mrrcv]
        findP(@view(result_mat[:, (i-1)*4+1]), max_element_each_row, @view(result_mat[:, (i-1)*4+3])) #, 1)
        result_mat[:, i*4] = ones(dim) * quantile(max_element_each_row, settings.alpha)
        i += 1
    end

    return (result_mat)

end


function DGP_boot(y, trend::Bool, bR, bU, eR, eU, lag, From_var, To_var, nb, M)

    T, num_y = size(y)

    # random number generator for drawing residuals
    Random.seed!(nb)
    rN = rand(1:T-lag, M - lag)  #a random vector of (M-lag) elements with each between 1 and T-lag
    yb = zeros(M, num_y)

    yb[1:lag, :] = @view(y[1:lag, :])               ## the first lag obs are the same as in the original data

    if trend==false
        num_x = 0
        xstar_U = [ones(M - lag) zeros(M - lag, num_y * lag)]
    else
        num_x = 1
        xstar_U = [ones(M - lag) @view(x[lag+1:M]) zeros(M - lag, num_y * lag)]
    end

    @inbounds for j ∈ lag+1:M  # j is the row index of yb

        @inbounds for i ∈ 1:lag
            xstar_U[j-lag, 1+num_x+(i-1)*num_y+1:1+num_x+i*num_y] = @view(yb[j-i, :])
        end

        @inbounds for endo ∈ 1:num_y
            if endo == To_var
                yb[j, endo] = dot(xstar_U[j-lag, :], bR[:, 1]) + eR[rN[j-lag]]
            else
                yb[j, endo] = dot(xstar_U[j-lag, :], bU[:, endo]) + eU[rN[j-lag], endo]
            end
        end
    end


    return yb
end

function findP(stat, mcv, result)  #, flag)

    dim = length(stat)
    L = length(mcv)
    #pvalue = Array{Float64}(undef, dim)

    @inbounds for i ∈ 1:dim
        ind = (mcv .> stat[i])
        result[i] = sum(ind) / L
    end


end

function var_soc(y, max_lag::Int64, deterministic::Symbol, IC="bic")

    T, num_y = size(y)
    nob_act = T - max_lag
    num_x = ifelse(deterministic == :trend, 1, 0)

    reg0 = _common_TS.get_reg(y, deterministic, max_lag)

    sq_num_y = num_y * num_y
    y_x = num_y * num_x

    y_reg = @view(y[max_lag+1:T, :])

    min_ic_value::Float64 = typemax(Float64)
    tbr::Int64 = 0

    @inbounds for j ∈ 1:max_lag

        reg = @view(reg0[max_lag+1:T, 1:1+num_x+num_y*j])

        b, e = var_ols(y_reg, reg)
        sigma = e' * e / size(e, 1)

        item1 = nob_act * log(det(sigma))
        item_2 = (sq_num_y * j + y_x)
        if IC == "aic"
            temp = item1 + 2 * item_2
        elseif IC == "hq"
            temp = item1 + 2 * log(log(nob_act)) * item_2
        else
            temp = item1 + log(nob_act) * item_2  #bic
        end
        if temp < min_ic_value
            min_ic_value = temp
            tbr = j
        end
    end
    return tbr
end

function var_ols(y, reg, fast=true)
    if fast == true
        b = (reg' * reg) \ (reg' * y)
    else
        b = reg \ y
    end
    e = y - reg * b
    return b, e
end

# function get_reg(y, x, lag, lag_first=true)
#     obs, num_y = size(y)  #bug: must be a matrix
#     #return的矩阵，实际从lag+1开始

#     if x != [0]
#         num_x = 1
#         reg = [ones(obs) @view(x[1:obs]) Array{Float64}(undef, obs, num_y * lag)]
#     else
#         num_x = 0
#         reg = [ones(obs) Array{Float64}(undef, obs, num_y * lag)]
#     end

#     if lag_first
#         @inbounds @simd for j ∈ 1:lag
#             reg[j+1:obs, 1+num_x+num_y*(j-1)+1:1+num_x+num_y*j] = @view(y[1:obs-j, :])
#         end
#     else
#         for s in 1:num_y
#             for j in 1:lag
#                 reg[j+1:obs, 1+num_x+num_y*(s-1)+j] = @view(y[1:obs-j, s])
#             end
#         end
#     end

#     return reg
# end



function var_ols_rest_b(y, lag, From_var, To_var, reg_full, trend::Bool)

    T, num_y = size(y)

    num_x = ifelse(trend, 1, 0)

    xstar_U = @view(reg_full[lag+1:T, 1:1+num_x+lag*num_y])

    ind = [1+num_x+From_var:num_y:1+num_x+(num_y-1)*lag+From_var;]
    L = size(xstar_U, 2) - lag  #L is the number of columns in xstar_R

    ystar_U = view(y, lag+1:T, :)
    ystar_R = view(y, lag+1:T, To_var)
    xstar_R = @view(xstar_U[:, setdiff(1:size(xstar_U, 2), ind)])

    bR = zeros(size(xstar_U, 2))
    rR = zeros(size(xstar_U, 1))
    bR0 = @view(bR[setdiff(1:end, ind)])

    bR0, rR = var_ols(ystar_R, xstar_R)



    bU = zeros(size(xstar_U, 2), num_y)
    rU = zeros(size(xstar_U, 1), num_y)

    for k ∈ 1:num_y
        if k != To_var
            # bU[:, k] = (xstar_U' * xstar_U) \ (xstar_U' * @view(y[lag+1:T, k:k]))
            # rU[:, k] = @view(y[lag+1:T, k:k]) - xstar_U * @view(bU[:, k])
            bU[:, k], rU[:, k] = var_ols(@view(y[lag+1:T, k:k]), xstar_U)

        end
    end

    return bR, bU, rR, rU

end

function bootstrap_RGC_MW(y, lag, From_var, To_var, reg_full, M, settings::testSettings)

    bR, bU, eR, eU = var_ols_rest_b(y, lag, From_var, To_var, reg_full, settings.trend)

    dim = M - settings.window_size + 1      #% T-swindow0+1   # in the main proc, M = swindow0 + 12 * 1 - 1 so, dim =12

    nboot = settings.n_boot

    TS_f = zeros(nboot, dim)
    TS_r = zeros(nboot, dim)
    TS_rr = zeros(nboot, dim)

    @inbounds for nb ∈ 1:nboot
        yb = DGP_boot(y, settings.trend , bR, bU, eR, eU, lag, From_var, To_var, nb, M)
        if settings.trend
            reg_full =_common_TS.get_reg(yb, :trend, lag + settings.d)
        else
            reg_full = _common_TS.get_reg(yb, :constant, lag + settings.d)
        end

        @inbounds for r2 ∈ settings.window_size:M
            rwGC = zeros(r2 - settings.window_size + 1)
            for r1 ∈ 1:(r2-settings.window_size+1)
                temp_a, temp_b = granger_cause_Mwald(yb, r1, r2, lag, From_var, To_var,  reg_full, settings)
                rwGC[r1] = temp_a
                if r1 == 1
                    TS_f[nb, r2-settings.window_size+1] = temp_a
                end
                if r1 == r2 - settings.window_size + 1
                    TS_r[nb, r2-settings.window_size+1] = temp_a
                end
            end
            if r2 == settings.window_size
                TS_rr[nb, r2-settings.window_size+1] = rwGC[1]
            else
                TS_rr[nb, r2-settings.window_size+1] = maximum(rwGC)
            end
        end

    end

    return vec(findmax(TS_f, dims=2)[1]), vec(findmax(TS_r, dims=2)[1]), vec(findmax(TS_rr, dims=2)[1])
end

function calculate_cov(reg, e, num_y, nob, k, d, robust)

    if robust == false
        sigma = e' * e / size(e, 1)
        iQT = inv(reg' * reg)
        cov = kron(sigma, iQT)

    else  ## robst == true
        reg_width = size(reg, 2)
        Q_hat = (reg' * reg) / nob
        In = Matrix{Float64}(I, num_y, num_y) #eye(num_y)
        V_hat = kron(In, Q_hat)
        IV_hat = inv(V_hat)

        #mxi = Array{Float64}(undef, nob, num_y * (num_y * (k + d) + 1 + num_x))
        mxi = Array{Float64}(undef, nob, num_y * reg_width)
        @inbounds @simd for i in 1:num_y
            mxi[:, (i-1)*reg_width+1:i*reg_width] = @view(e[:, i]) .* @view(reg[:, 1:reg_width])
        end

        # for i ∈ 1:nob
        #     mxi[i, :] = kron(@view(e[i, :]), @view(reg[i, :]))
        # end
        W_hat = (mxi' * mxi) / nob
        cov = IV_hat * W_hat * IV_hat

    end

    return (cov)
end

function granger_cause_Mwald(y, start_t, end_t, lag, FROM, TO, reg_full, settings::testSettings)

    k = lag
    d = settings.d
    nob = (end_t - start_t + 1) - k - d
    T, num_y = size(y)


    num_x = ifelse(settings.trend == true, 1, 0)

    reg = @view(reg_full[start_t+k+d:end_t, :])


    ystar = @view(y[start_t+k+d:end_t, :])
    #writedlm( "reg.csv",  reg, ',')

    b, e = var_ols(ystar, reg)
    cov = calculate_cov(reg, e, num_y, nob, k, d, settings.robust)
    se_b = sqrt.(diag(cov) / nob)

    d1, d2 = size(b)

    vb = reshape(b, d1 * d2, 1)


    R = zeros(k, d1 * d2)

    for j ∈ 1:k
        R[j, (TO-1)*d1+1+num_x+num_y*(j-1)+FROM] = 1
    end

    scoeff = sum(R * (vb ./ se_b))

    if settings.robust
        wald = nob * (R * vb)' * ((R * cov * R') \ (R * vb))
    else
        wald = (R * vb)' * ((R * cov * R') \ (R * vb))
    end






    return wald[1], scoeff
end



end