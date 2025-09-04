using CSV, DataFrames 
using DelimitedFiles
using Plots
using Printf
using LaTeXStrings

N = 10000
q = 3
alpha = 0.1
c0 = 1
he = 0.
hp = 0.01
yearvariant = 2025

filename = "EPO_$yearvariant _N$(N)_q$(q)_alpha$(@sprintf("%.2f", alpha))_c0$(@sprintf("%.2f", c0))_he$(@sprintf("%.2f", he))_hp$(@sprintf("%.2f", hp)).txt"

#df = CSV.read(joinpath("qv_EPO", filename), DataFrame; delim='\t')  
df = readdlm(joinpath("qv_EPO_$yearvariant", filename), '\t')  # For CSV-style files
# print(df)

p = df[:, 1]
cE = df[:, 2]
cP = df[:, 3]
diss = df[:, 4]


# plot(p, cE, seriestype=:scatter)
# xlims!(0, 1)
# ylims!(0, 1)
# xlabel!("\$p\$")
# ylabel!("\$c_S\$")
# title!("alpha = $(@sprintf("%.2f", alpha)); he = $(@sprintf("%.2f", he)); hp = $(@sprintf("%.2f", hp))")


# plot(p, cP, seriestype=:scatter)
# xlims!(0, 1)
# ylims!(0, 1)
# xlabel!("\$p\$")
# ylabel!("\$c_sigma\$")
# title!("alpha = $(@sprintf("%.2f", alpha)); he = $(@sprintf("%.2f", he)); hp = $(@sprintf("%.2f", hp))")

# plot(p, diss, seriestype=:scatter)
# xlims!(0, 1)
# ylims!(0, 1)
# xlabel!("\$p\$")
# ylabel!("\$d\$")
# title!("alpha = $(@sprintf("%.2f", alpha)); he = $(@sprintf("%.2f", he)); hp = $(@sprintf("%.2f", hp))")
marksize = 4

plot()
plot!(p, cE, seriestype=:scatter, shape =:circle, ms = marksize, label = "Expressed")
plot!(p, cP, seriestype=:scatter, shape = :utriangle, ms = marksize, label = "Private")
xlims!(0, 1)
ylims!(0, 1)
xlabel!("\$p\$")
ylabel!("\$c\$")
title!("alpha = $(@sprintf("%.2f", alpha)); he = $(@sprintf("%.2f", he)); hp = $(@sprintf("%.2f", hp))")