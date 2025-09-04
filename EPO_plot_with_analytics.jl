using CSV, DataFrames 
using DelimitedFiles
using Plots
using Plots.Measures
using Printf
using LaTeXStrings

N = 10000
q = 3
alpha = 0.1
c0 = 1.
yearvariant = "2018"

plot_name = "EPO_$yearvariant _N$(N)_q$(q)_alpha$(@sprintf("%.2f", alpha))_c0$(@sprintf("%.2f", c0))_with_MFA.png"

# h = [[0, 0], [0.01, 0], [0.05, 0]]
# h = [[0, 0], [0.01, 0], [0.05, 0], [0, 0.01], [0, 0.05]]
h = [[0.01, 0.0]]

p1 = plot(xlabel="\$p\$", ylabel="\$c_S\$", xlims=(0, 1), ylims=(0, 1), margin=5mm)
p2 = plot(xlabel="\$p\$", ylabel="\$c_σ\$", xlims=(0, 1), ylims=(0, 1), margin=5mm)
p3 = plot(xlabel="\$p\$", ylabel="\$d\$", xlims=(0, 1), ylims=(0, 1), margin=5mm)

for htemp in h
    he, hp = htemp
    filename = "EPO_$yearvariant _N$(N)_q$(q)_alpha$(@sprintf("%.2f", alpha))_c0$(@sprintf("%.2f", c0))_he$(@sprintf("%.2f", he))_hp$(@sprintf("%.2f", hp)).txt"

    #df = CSV.read(joinpath("qv_EPO", filename), DataFrame; delim='\t')  
    df = readdlm(joinpath("qv_EPO_$yearvariant", filename), '\t')  # For CSV-style files
    # print(df)

    p = df[:, 1]
    cE = df[:, 2]
    cP = df[:, 3]
    diss = df[:, 4]

    pathanalytics = "C:/Users/Basia/OneDrive - Politechnika Wroclawska/Pulpit/PhD/Current work/ExternalInfluence/MFA"
    # pathanalytics = "qv_EPO_$yearvariant"
    filenameanalytics = "EPO_$(yearvariant)_q$(q)_alpha_$(@sprintf("%.2f", alpha))_h_e_$(@sprintf("%.2f", he))_h_p_$(@sprintf("%.2f", hp)).txt"

    dfanalytics = readdlm(joinpath(pathanalytics, filenameanalytics), '\t')
    panalytics = dfanalytics[:, 1]
    cEanalytics = dfanalytics[:, 2]
    cPanalytics = dfanalytics[:, 3]
    dissanalytics = dfanalytics[:, 4]


    msize = 2
    plot!(p1, p, cE, seriestype=:scatter, markersize = msize, framestyle=:box, xlabelfontsize=14, ylabelfontsize=14,
     label = "he = $(@sprintf("%.2f", he)); hp = $(@sprintf("%.2f", hp))")
    plot!(p2, p, cP, seriestype=:scatter, markersize = msize, framestyle=:box, xlabelfontsize=14, ylabelfontsize=14, label = false)
    plot!(p3, p, diss, seriestype=:scatter, markersize = msize, framestyle=:box, xlabelfontsize=14, ylabelfontsize=14, label = false)

    plot!(p1, panalytics, cEanalytics,framestyle=:box, xlabelfontsize=14, ylabelfontsize=14,
     label = "he = $(@sprintf("%.2f", he)); hp = $(@sprintf("%.2f", hp))")
    plot!(p2, panalytics, cPanalytics, framestyle=:box, xlabelfontsize=14, ylabelfontsize=14, label = false)
    plot!(p3, panalytics, dissanalytics, framestyle=:box, xlabelfontsize=14, ylabelfontsize=14, label = false)


end

# xlims!(0, 1)
# ylims!(0, 1)
# xlabel!("\$p\$")
# ylabel!("\$c_S\$")
# title!("model_$yearvariant, alpha = $(@sprintf("%.2f", alpha))")


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

final_plot = plot(p1, p2, p3, layout=(1, 3), size=(1200, 400), plot_title="Model $yearvariant, q=$q, alpha = $(@sprintf("%.2f", alpha))", dpi=300)
plot!()
# savefig(joinpath("qv_EPO_figures/", plot_name))
