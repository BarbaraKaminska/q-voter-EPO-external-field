using CSV, DataFrames 
using DelimitedFiles
using Plots
using Printf
using LaTeXStrings
using Plots: mm

N = 1000000
q = 3
p = 0.25
alpha = 0.1
# c0 = 0.5
c0_ =  0.:0.1:1.

yearvariant = "2025"

println("Start plot")


colors = [:black, :green, :blue]
shapes = [:circle, :diamond, :utriangle]

# plot_name = "EPO_traj_$yearvariant _N$(N)_q$(q)_p$(@sprintf("%.2f", p))_alpha$(@sprintf("%.2f", alpha))_c0$(@sprintf("%.2f", c0)).png"

# h = [[0, 0], [0.01, 0], [0.05, 0], [0, 0.01], [0, 0.05]]
h = [[0., 0.], [0.0, 0.05], [0.05, 0.0]]
plot(dpi=300)

for (idx_c, c0) in enumerate(c0_)
    # println(c0)

    filename = "EPO_traj_$yearvariant _N$(N)_q$(q)_p$(@sprintf("%.2f", p))_alpha$(@sprintf("%.2f", alpha))_c0$(@sprintf("%.2f", c0)).txt"
    # filename = "EPO_traj_2025 _N1000000_q3_p0.20_alpha0.10_c00.40.txt"
    #df = CSV.read(joinpath("qv_EPO", filename), DataFrame; delim='\t')  
    df = readdlm(joinpath("qv_EPO_$yearvariant", filename), '\t')  # For CSV-style files
    # print(df)

    t = df[:, 1]
    he = df[:, 2]
    hp = df[:, 3]        
    cE = df[:, 4]
    cP = df[:, 5]
    diss = df[:, 6]

    global tmax = t[end]
    
    if (idx_c == 1)
        plot!(t, cE, seriestype=:line, color = colors[1], label = "Expressed")
        plot!(t, he, seriestype=:line, color = colors[2], label = "External field expressed")
        plot!(t, hp, seriestype=:line, color = colors[3], label = "External field private", legend=:bottomright)
    else
        plot!(t, cE, seriestype=:line, color = colors[1], label = "")
        plot!(t, he, seriestype=:line, color = colors[2], label = "")
        plot!(t, hp, seriestype=:line, color = colors[3], label = "")
        # plot!(t, cE, seriestype=:scatter, color = colors[idx], markershape = shapes[idx],
        # label="", size=(650, 400), left_margin = 2Plots.mm, right_margin=2Plots.mm, legend=false)
    end
    
end

println("Plot ready")
xlims!(0, tmax+10)
ylims!(0, 1.)
xlabel!("time")
ylabel!("Fraction of agents \n with ecofriendly behavior")

# title!("model_$yearvariant, alpha = $(@sprintf("%.2f", alpha))")
title!("")

plotname = "trajectory $yearvariant _N$(N)_q$(q)_p$(@sprintf("%.2f", p))_alpha$(@sprintf("%.2f", alpha)).png"


savefig(joinpath("qv_EPO_figures_trajectories/", plotname))
plot!()
