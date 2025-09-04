include("../EPO_functions.jl")

yearvariant = "2025"

N = 1000000
q = 3
p = 0.2
alpha = 0.1
# c0 = .5
# he = 0.05
# hp = 0.
MCS_term = 10000
MCS_total = 500
MCS_step = 1
# h = [[0, 0], [0.01, 0], [0.05, 0], [0, 0.01], [0, 0.05]]
h = [[0, 0], [0., 0.01], [0.01, 0]]
c0_ = 0.:0.1:1.
println("Start")



for c0 in c0_

    for htemp in h
        he, hp = htemp
   
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
            if mcs % MCS_step == 0 
                exporttrajectory(N, q, alpha, c0, he, hp, p, mcs, cE, cP, diss, yearvariant)
            end
            expressedopinions, privateopinions = singleupdate_2025!(expressedopinions, privateopinions, N, q, p, alpha, he, hp)
        end
    end
    # println("p = $(@sprintf("%.2f", p));\tExpressed = $(@sprintf("%.2f", cE));\t Private = $(@sprintf("%.2f", cP));\t Dissonance = $(@sprintf("%.2f", diss))")

end

println("Done")