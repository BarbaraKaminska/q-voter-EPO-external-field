using CSV, DataFrames 
using DelimitedFiles
using Plots
using Printf
using LaTeXStrings

N = 10000
q = 3
alpha = 0.1
c0 = 1.
yearvariant = "2025"

plot_name = "EPO_$yearvariant _N$(N)_q$(q)_alpha$(@sprintf("%.2f", alpha))_c0$(@sprintf("%.2f", c0)).png"

# h = [[0, 0], [0.01, 0], [0.05, 0], [0, 0.01], [0, 0.05]]
h = [[0, 0], [0.01, 0], [0, 0.01]]
# h = [[0.0, 0.01]]
plot()

for htemp in h
    he, hp = htemp
    filename = "EPO_$(yearvariant)_N$(N)_q$(q)_alpha$(@sprintf("%.2f", alpha))_c0$(@sprintf("%.2f", c0))_he$(@sprintf("%.2f", he))_hp$(@sprintf("%.2f", hp)).txt"

    #df = CSV.read(joinpath("qv_EPO", filename), DataFrame; delim='\t')  
    df = readdlm(joinpath("qv_EPO_$yearvariant", filename), '\t')  # For CSV-style files
    # print(df)

    p = df[:, 1]
    cE = df[:, 2]
    cP = df[:, 3]
    diss = df[:, 4]

    if he == 0.0 && hp == 0.0
        plot!(p, cE, seriestype=:scatter, markershape = :circle, msize = 4, label = "No external field")
    elseif he != 0.0 && hp == 0.0
        plot!(p, cE, seriestype=:scatter, markershape = :diamond, msize = 4, label = "he = $(@sprintf("%.2f", he)); hp = $(@sprintf("%.2f", hp))")
    else 
        plot!(p, cE, seriestype=:scatter, markershape = :square, msize = 3, label = "he = $(@sprintf("%.2f", he)); hp = $(@sprintf("%.2f", hp))")
    end


end

xlims!(0, 1)
ylims!(0, 1.05)
xlabel!("\$p\$")
ylabel!("\$c_S\$")
title!("model_$yearvariant, alpha = $(@sprintf("%.2f", alpha))")
# c0 = 0.5

# for htemp in h
#     he, hp = htemp
#     filename = "EPO_$yearvariant _N$(N)_q$(q)_alpha$(@sprintf("%.2f", alpha))_c0$(@sprintf("%.2f", c0))_he$(@sprintf("%.2f", he))_hp$(@sprintf("%.2f", hp)).txt"

#     #df = CSV.read(joinpath("qv_EPO", filename), DataFrame; delim='\t')  
#     df = readdlm(joinpath("qv_EPO_$yearvariant", filename), '\t')  # For CSV-style files
#     # print(df)

#     p = df[:, 1]
#     cE = map((x) -> 0.5 + abs(0.5-x),  df[:, 2])
#     cP = map((x) -> 0.5 + abs(0.5-x),  df[:, 3])
#     diss = df[:, 4]


#     plot!(p, cE, seriestype=:scatter, label = nothing)
#     xlims!(0, 1)
#     ylims!(0, 1)
#     xlabel!("\$p\$")
#     ylabel!("\$c_S\$")
#     title!("alpha = $(@sprintf("%.2f", alpha))")

# end

plot!(dpi=300)
# savefig(joinpath("qv_EPO_figures_$yearvariant/", plot_name))
