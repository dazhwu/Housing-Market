#using MKL
using Dates
using DataFrames
using Plots
using CSV
using HypothesisTests

include("./lib/tvgc.jl")
include("./lib/exuber.jl")
include("./lib/exuber_old.jl")
include("./lib/common_proc.jl")


df0 = CSV.read("bubble.csv", DataFrame)
bsadf, cv, d_label=gsadf.exuber(df0,1,0)

col_width=size(df0,2)
for i in 2:col_width

    # wboot=exuber_old.wmboot(ts,swindow0,2,adflag,Tb);
    # Q_boot=wboot[2] * ones(length(bsadf),1);
    
    #Q_boot=wboot[2] * ones(length(results.bsadf),1);
    
    #
    #Tb = 4 * 2 + swindow0 - 1
    #s_b=gsadf.wmboot(ts,swindow0, adflag,  499)
 
    
    plot(d_label, bsadf[:,i-1], label="bsadf", dpi=1200, lw=2, linecolor=:black)
    
    plot!(d_label, cv[:,i-1], label="95% cv", line=(:dot, 2), linecolor=:red, dpi=1200)
    
    #plot!(d_label, s_b[:, 2], label="wild cv", line=(:dot, 2), linecolor=:blue, dpi=1200)
    
    savefig("./png/Bubbles/" * names(df0)[i] * ".png")
    # plot(d_label_old, bsadf, label="bsadf", dpi=1200, lw=2, linecolor=:black)
    # plot!(d_label_old, Q_boot, label="95% wmboot cv", line=(:dot, 2), linecolor=:blue, dpi=1200)
    # savefig("./png/Bubbles_old/" * MSA_name * ".png")

end


