using Random
using Printf
using StatsBase
using DelimitedFiles

function c(X)
    return 0.5*(mean(X)+1)
end

function d(X, Y)
    return 0.5*(1 - mean(X .* Y))
end

function exportresults(N, q, alpha, c0, he, hp, p, expressed, private, dissonance, yearvariant)
    filename = "EPO_$yearvariant _N$(N)_q$(q)_alpha$(@sprintf("%.2f", alpha))_c0$(@sprintf("%.2f", c0))_he$(@sprintf("%.2f", he))_hp$(@sprintf("%.2f", hp)).txt"
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
    filename = "EPO_traj_$yearvariant _N$(N)_q$(q)_p$(@sprintf("%.2f", p))_alpha$(@sprintf("%.2f", alpha))_c0$(@sprintf("%.2f", c0))_he$(@sprintf("%.2f", he))_hp$(@sprintf("%.2f", hp)).txt"
    # println(filename)
    open(joinpath("qv_EPO_$yearvariant", filename), "a") do io
        writedlm(io,  [t expressed private dissonance], '\t')
    end
end

function exporttrajectorywithfield(N, q, alpha, c0, he, hp, p, t, expressed, private, dissonance, yearvariant)
    filename = "EPO_traj_$yearvariant _N$(N)_q$(q)_p$(@sprintf("%.2f", p))_alpha$(@sprintf("%.2f", alpha))_c0$(@sprintf("%.2f", c0)).txt"
    # println(filename)
    open(joinpath("qv_EPO_$yearvariant", filename), "a") do io
        writedlm(io,  [t he hp expressed private dissonance], '\t')
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


# Model with avoiding dissonance
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

# model withouth avoiding dissonance
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


