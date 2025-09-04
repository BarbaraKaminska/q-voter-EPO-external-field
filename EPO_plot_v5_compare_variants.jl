using CSV, DataFrames 
using DelimitedFiles
using Plots
using Plots.Measures
using Printf
using LaTeXStrings

N = 10000
q = 3
alpha = 0.1
c0 = 1

h = [0.05, 0.0]
he, hp = h

plot_name = "EPO_N$(N)_q$(q)_alpha$(@sprintf("%.2f", alpha))_c0$(@sprintf("%.2f", c0))_he$(@sprintf("%.2f", he))_hp$(@sprintf("%.2f", hp)).png"


yearvariants = ["2024_v0", "2024_v1", "2024_v2"]

p1 = plot(xlabel="\$p\$", ylabel="\$c_S\$", xlims=(0, 1), ylims=(0, 1), margin=5mm)
p2 = plot(xlabel="\$p\$", ylabel="\$c_σ\$", xlims=(0, 1), ylims=(0, 1), margin=5mm)
p3 = plot(xlabel="\$p\$", ylabel="\$d\$", xlims=(0, 1), ylims=(0, 1), margin=5mm)

for yearvariant in yearvariants
    
    filename = "EPO_$yearvariant _N$(N)_q$(q)_alpha$(@sprintf("%.2f", alpha))_c0$(@sprintf("%.2f", c0))_he$(@sprintf("%.2f", he))_hp$(@sprintf("%.2f", hp)).txt"

    #df = CSV.read(joinpath("qv_EPO", filename), DataFrame; delim='\t')  
    df = readdlm(joinpath("qv_EPO_$yearvariant", filename), '\t')  # For CSV-style files
    # print(df)

    p = df[:, 1]
    cE = df[:, 2]
    cP = df[:, 3]
    diss = df[:, 4]


    plot!(p1, p, cE, seriestype=:scatter, framestyle=:box, xlabelfontsize=14, ylabelfontsize=14,
     label = "$(yearvariant[end-1:end])")
    plot!(p2, p, cP, seriestype=:scatter, framestyle=:box, xlabelfontsize=14, ylabelfontsize=14, label = false)
    plot!(p3, p, diss, seriestype=:scatter, framestyle=:box, xlabelfontsize=14, ylabelfontsize=14, label = false)

end



final_plot = plot(p1, p2, p3, layout=(1, 3), size=(1200, 400), plot_title="Model $yearvariant, q=$q, alpha = $(@sprintf("%.2f", alpha))", dpi=300)
plot!()
savefig(joinpath("qv_EPO_figures/", plot_name))
