include("../EPO_functions.jl")

yearvariant = "2025"

N = 1000000
q = 3
p = 0.25
alpha = 0.1
# c0 = .5
he = 0.02
hp = 0.02
MCS_he = 2000
MCS_hp = 2000
MCS_none = 2000



c0_ = 0.:0.1:1.
println("Start")

he_, hp_ = externalfield(he, hp, MCS_he, MCS_hp, MCS_none) 
he_ = vcat(he_, externalfield(he, hp, MCS_he, MCS_hp, MCS_none)[1])
hp_ = vcat(hp_, externalfield(he, hp, MCS_he, MCS_hp, MCS_none)[2])

MCS_total = length(he_) - 1 
MCS_step = 1

for c0 in c0_

    global cE = 0
    global cP = 0
    global diss = 0
    expressedopinions = ones(N)

    for i in 1:1:Int(round(N*(1-c0)))
        expressedopinions[i] = (-1)
    end
    
    privateopinions = copy(expressedopinions)

    for mcs=0:1:MCS_total
        # println(mcs)
        cE = c(expressedopinions)
        cP = c(privateopinions)
        diss = d(expressedopinions, privateopinions)

        he_temp = he_[mcs + 1]
        hp_temp = hp_[mcs + 1]        
        if mcs % MCS_step == 0 
            exporttrajectorywithfield(N, q, alpha, c0, he_temp, hp_temp, p, mcs, cE, cP, diss, yearvariant)
        end
        expressedopinions, privateopinions = singleupdate_2025!(expressedopinions, privateopinions, N, q, p, alpha, he_temp, hp_temp)
    end

    # println("p = $(@sprintf("%.2f", p));\tExpressed = $(@sprintf("%.2f", cE));\t Private = $(@sprintf("%.2f", cP));\t Dissonance = $(@sprintf("%.2f", diss))")

end

println("Done")