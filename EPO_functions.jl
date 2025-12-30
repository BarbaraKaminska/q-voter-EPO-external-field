using Random
using Printf
using StatsBase
using DelimitedFiles
using Symbolics
using NLsolve
using LinearAlgebra


function c(X)
    return 0.5*(mean(X)+1)
end

function d(X, Y)
    return 0.5*(1 - mean(X .* Y))
end

function exportresults(N, q, alpha, c0, he, hp, p, expressed, private, dissonance, yearvariant)
    filename = "EPO_$(yearvariant)_N$(N)_q$(q)_alpha$(@sprintf("%.2f", alpha))_c0$(@sprintf("%.2f", c0))_he$(@sprintf("%.2f", he))_hp$(@sprintf("%.2f", hp)).txt"
    # println(filename)
    open(joinpath("qv_EPO_$yearvariant", filename), "a") do io
        writedlm(io,  [p expressed private dissonance], '\t')
        # for val in [p, expressed, private, dissonance]
        #     @printf(io, "%.4f ", val)
        # end
        # @printf(io, "\n")
    end
end

function exporttrajectory(N, q, alpha, c0, he, hp, p, t, expressed, private, dissonance, yearvariant)
    filename = "EPO_traj_$(yearvariant)_N$(N)_q$(q)_p$(@sprintf("%.2f", p))_alpha$(@sprintf("%.2f", alpha))_c0$(@sprintf("%.2f", c0))_he$(@sprintf("%.2f", he))_hp$(@sprintf("%.2f", hp)).txt"
    # println(filename)
    open(joinpath("qv_EPO_$yearvariant", filename), "a") do io
        writedlm(io,  [t expressed private dissonance], '\t')
    end
end

function exporttrajectorywithfield(N, q, alpha, c0, he, hp, p, t, expressed, private, dissonance, yearvariant)
    filename = "EPO_traj_$(yearvariant)_N$(N)_q$(q)_p$(@sprintf("%.2f", p))_alpha$(@sprintf("%.2f", alpha))_c0$(@sprintf("%.2f", c0)).txt"
    # println(filename)
    open(joinpath("qv_EPO_$yearvariant", filename), "a") do io
        writedlm(io,  [t he hp expressed private dissonance], '\t')
    end
end

function exportwithstability(filename, df, yearvariant)
    println(filename)
    open(joinpath("qv_EPO_$yearvariant", filename), "a") do io
        writedlm(io,  [df[:, 1] df[:, 2] df[:, 3] df[:, 4] df[:, 5] df[:, 6]], '\t')
    end
end

# Generate variable external field
function externalfield(he, hp, t_he, t_hp, t_none)
    he_ = zeros(t_none + t_he + t_none + t_hp + 1)
    hp_ = zeros(t_none + t_he + t_none + t_hp + 1)
    he_[t_none+1:t_none+t_he] .= he
    hp_[2*t_none+t_he+1:end] .= hp
    return he_, hp_
end


# Model without self-anticonformity
function singleupdate_2025!(expr, priv, N, q, p, alpha, he, hp)
    for _ = 1:1:N
        target = rand(1:N)
        if rand() < alpha
            # Update private opinion
            if rand() < p
                # independence
                if rand() < 0.5
                    priv[target] *=(-1)
                end
            else
                # conformity
                if rand() < hp
                    priv[target] = 1
                else
                    q_panel = sample(1:N, q; replace = false)
                    if sum(expr[q_panel]) == q*expr[target]
                        priv[target] = expr[q_panel[1]]
                    end
                end
            end
        else        
            # Update expressed opinion
            if rand() < p
                # independence
                expr[target] = priv[target]
            else 
                # conformity
                if rand() < he
                    expr[target] = 1
                else
                    q_panel = sample(1:N, q; replace = false)
                    if expr[target] == priv[target]
                        # compliance 
                        if abs(sum(expr[q_panel])) == q 
                            expr[target] = expr[q_panel[1]]
                        end
                    else 
                        for i in q_panel
                            if expr[i] == priv[target]
                                expr[target] = priv[target]
                                break
                            end
                        end
                    end
                end
            end 
        end
    end
    return expr, priv
end

# model with self-anticonformity
function singleupdate_2018!(expr, priv, N, q, p, alpha, he, hp)
    for _ = 1:1:N
        target = rand(1:N)
        if rand() < alpha
            # Update private opinion
            if rand() < p
                # independence
                if rand() < 0.5
                    priv[target] *=(-1)
                end
            else
                # conformity
                if rand() < hp
                    priv[target] = 1
                else
                    q_panel = sample(1:N, q; replace = false)
                    if abs(sum(expr[q_panel])) == q 
                        priv[target] = expr[q_panel[1]]
                    end
                end
            end
        else        
            # Update expressed opinion
            if rand() < p
                # independence
                expr[target] = priv[target]
            else 
                # conformity
                if rand() < he
                    expr[target] = 1
                else
                    q_panel = sample(1:N, q; replace = false)
                    if expr[target] == priv[target]
                        # compliance 
                        if abs(sum(expr[q_panel])) == q 
                            expr[target] = expr[q_panel[1]]
                        end
                    else 
                        for i in q_panel
                            if expr[i] == priv[target]
                                expr[target] = priv[target]
                                break
                            end
                        end
                    end
                end
            end 
        end
    end
    return expr, priv
end

# stability check - model with self-anticonformity
function jacobian_2025()
    @variables uu ud du dd α p he hp q

    # Define the symbolic system equations
    F1 = uu*(1 - α*p/2 - (1-α)*(1-p)*(1-he)*(1-uu-ud)^q) +
        ud*α*(p/2 + (1-p)*hp + (1-p)*(1-hp)*(uu+ud)^q) +
        du*(1-α)*(1 - (1-p)*(1-he)*(1-uu-ud)^q) #- uu

    F2 = uu*α*p/2 +
        ud*((1-α)*(1-p)*(he + (1-he)*(uu+ud)^q) + α*(1-p/2) - α*(1-p)*(hp + (1-hp)*(uu+ud)^q)) +
        dd*(1-α)*(1-p)*(he + (1-he)*(uu+ud)^q) #- ud

    F3 = uu*(1-α)*(1-p)*(1-he)*(1-uu-ud)^q +
        du*((1-α)*(1-p)*(1-he)*(1-uu-ud)^q + α*(1 - p/2 - (1-p)*hp - (1-p)*(1-hp)*(1-uu-ud)^q)) +
        dd*α*(p/2 + (1-p)*hp) # - du

    F4 = 1 - (F1 + F2 + F3)

    # Create the system vector and variable vector
    F_vec = [F1, F2, F3, F4]
    vars = [uu, ud, du, dd]

    # Calculate the symbolic Jacobian
    # J_symbolic = Symbolics.jacobian(F_vec, vars)


    J = Symbolics.jacobian(F_vec, vars)

    # === Optional: simplify the Jacobian ===
    J_simplified = simplify.(J)
    return J_simplified

end


function jacobian_2018()
    @variables uu ud du dd α p he hp q

    # Define the symbolic system equations
    F1 = uu*(1 - α*p/2 - α*(1-p)*(1-hp)*(1-uu-ud)^q - (1-α)*(1-p)*(1-he)*(1-uu-ud)^q) +
        ud*α*(p/2 + (1-p)*hp + (1-p)*(1-hp)*(uu+ud)^q) +
        du*(1-α)*(1 - (1-p)*(1-he)*(1-uu-ud)^q) #- uu

    F2 = uu*(α*p/2 +α*(1-p)*(1-hp)*(1-uu-ud)^q ) +
        ud*((1-α)*(1-p)*(he + (1-he)*(uu+ud)^q) + α*(1-p/2) - α*(1-p)*(hp + (1-hp)*(uu+ud)^q)) +
        dd*(1-α)*(1-p)*(he + (1-he)*(uu+ud)^q) #- ud

    F3 = uu*(1-α)*(1-p)*(1-he)*(1-uu-ud)^q +
        du*((1-α)*(1-p)*(1-he)*(1-uu-ud)^q + α*(1 - p/2 - (1-p)*hp - (1-p)*(1-hp)*(1-uu-ud)^q)) +
        dd*α*(p/2 + (1-p)*hp + (1-p)*(1-hp)*(uu+ud)^q) # - du

    F4 = 1 - (F1 + F2 + F3)

    # Create the system vector and variable vector
    F_vec = [F1, F2, F3, F4]
    vars = [uu, ud, du, dd]

    # Calculate the symbolic Jacobian
    # J_symbolic = Symbolics.jacobian(F_vec, vars)


    J = Symbolics.jacobian(F_vec, vars)

    # === Optional: simplify the Jacobian ===
    J_simplified = simplify.(J)
    return J_simplified

end

function checkstability(alpha_, p_, he_, hp_, q_, uu_, ud_, du_, dd_, jacobian)
    @variables uu ud du dd α p he hp q

    # Substitute fixed point + parameters into symbolic Jacobian
    subs_dict = Dict(α => alpha_, p => p_, he => he_, hp => hp_, q => q_,
        uu => uu_, ud => ud_,
        du => du_, dd => dd_)

    numericjacobian = substitute(jacobian, subs_dict) |> eval

    # Calculate eigenvalues
    eigenvalues = eigvals(numericjacobian)
    # println("Eigenvalues: ", eigenvalues)

    # Stability analysis
    if all(abs.(eigenvalues) .< 1)
        # println("Fixed point is STABLE")
        return true
    else
        # println("Fixed point is UNSTABLE")
        return false
    end
    
end
