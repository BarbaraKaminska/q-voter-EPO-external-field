include("EPO_functions.jl")

yearvariant = "2018"
q_ = [3, 5]
alpha_ = [0.1, 0.5]

# h = [[0, 0]]
h = [[0, 0], [0.01, 0], [0.05, 0], [0, 0.01], [0, 0.05]]

if yearvariant == "2025"
    jacobian_func = jacobian_2025()
    println("Jacobian for 2025 done")
elseif yearvariant == "2018"
    jacobian_func = jacobian_2018()
    println("Jacobian for 2018 done")
else 
    error("Unknown year variant: $yearvariant")
end

println(pwd())

for q in q_
for alpha in alpha_
for htemp in h
    he, hp = htemp


    println("he = $(@sprintf("%.2f", he)); hp = $(@sprintf("%.2f", hp))")

    pathanalytics = "C:/Users/Basia/OneDrive - Politechnika Wroclawska/Pulpit/PhD/Current work/ExternalInfluence/MFA"
    # cd(pathanalytics)
    # print("Current directory: ", pwd())
    # foreach(readdir()) do f
    #     println("\nObject: ", f)
    #     # dump(stat(f)) # you can customize what you want to print
    # end
        file = "EPO_$(yearvariant)_q$(q)_alpha_$(@sprintf("%.2f", alpha))_h_e_$(@sprintf("%.2f", he))_h_p_$(@sprintf("%.2f", hp))_4vars.txt"
        try
        # Load each file into a DataFrame
        df_ = readdlm(joinpath(pathanalytics, file), '\t')
        # append a column of zeros
        df_ = hcat(df_, zeros(size(df_, 1)))
        # println(df_[1:5, :])

        for rowid in eachindex(df_[:, 1])
            p_temp = df_[rowid, 1]
            uu_temp = df_[rowid, 2]
            ud_temp = df_[rowid, 3]
            du_temp = df_[rowid, 4]
            dd_temp = df_[rowid, 5]

            if checkstability(alpha, p_temp, he, hp, q, uu_temp, ud_temp, du_temp, dd_temp, jacobian_func)
                df_[rowid, 6] = 1.0
            else
                df_[rowid, 6] = 0.0
            end

        end

        # println("Stability check done.")
        # save the updated dataframe
        filetosave = "EPO_$(yearvariant)_q$(q)_alpha_$(@sprintf("%.2f", alpha))_h_e_$(@sprintf("%.2f", he))_h_p_$(@sprintf("%.2f", hp))_with_stability.txt"
        exportwithstability(filetosave, df_, yearvariant)

        println("Done")

    catch e
        # println("File not found for he = $(@sprintf("%.2f", he)), hp = $(@sprintf("%.2f", hp))")
        println(e)
        continue
    end
end
end
end