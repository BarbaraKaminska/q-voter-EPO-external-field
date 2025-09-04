include("../EPO_functions.jl")

yearvariant = "2018"

N = 10000
q = 3
alpha = 0.1
c0 = 1.0
# he = 0.05
# hp = 0.
MCS_term = 10000
MCS_total = 12000
MCS_step = 100

println("Start")

# h = [[0, 0], [0.01, 0], [0.05, 0], [0, 0.01], [0, 0.05]]
h = [[0., 0.01], [0.01, 0]]

for htemp in h
    he, hp = htemp

    for p = 0:0.01:1
        global cE = 0
        global cP = 0
        global diss = 0
        nmeasure = 0
        expressedopinions = ones(N)

        for i in 1:1:N
            if rand() > c0
                expressedopinions[i] = (-1)
            end
        end
        
        privateopinions = copy(expressedopinions)

        for mcs=1:1:MCS_total
            expressedopinions, privateopinions = singleupdate_2018!(expressedopinions, privateopinions, N, q, p, alpha, he, hp)
            # println(mcs)
            if mcs >= MCS_term && mcs % MCS_step == 0
                cE += c(expressedopinions)
                cP += c(privateopinions)
                diss += d(expressedopinions, privateopinions)
                nmeasure += 1
            end
        end

        cE/=nmeasure
        cP/=nmeasure
        diss/=nmeasure
        println("p = $(@sprintf("%.2f", p));\tExpressed = $(@sprintf("%.2f", cE));\t Private = $(@sprintf("%.2f", cP));\t Dissonance = $(@sprintf("%.2f", diss))")
        exportresults(N, q, alpha, c0, he, hp, p, cE, cP, diss, yearvariant)
    end

end