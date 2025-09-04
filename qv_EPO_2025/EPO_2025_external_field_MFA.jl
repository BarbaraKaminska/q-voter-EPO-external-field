using NLsolve
using Printf

# --- Parameters ---
p_range = 0.0:0.001:1.0
alpha_range = [0.1, 0.5]       # or 0.1:0.1:0.9
q_range = [3]
h_p_range = [0.01]
h_e_range = [0.]

# --- Main loops ---
for h_p in h_p_range
    for h_e in h_e_range
        for q in q_range
            for alpha in alpha_range

                cSUp = Float64[]
                cSigmaUp = Float64[]
                dissonance = Float64[]

                p_h = Float64[]
                cSup_h = Float64[]
                cSigmaUp_h = Float64[]
                dissonance_h = Float64[]

                # Output file
                baseFileName = @sprintf(
                    "EPO_2025_q%d_alpha_%.2f_h_e_%.2f_h_p_%.2f.txt",
                    q, alpha, h_e, h_p
                )
                export_file = joinpath("qv_EPO_2025", baseFileName)

                for p in p_range
                    # System of equations
                    function system!(F, x)
                        uu, ud, du, dd = x
                        F[1] = uu*(1 - alpha*p/2 - (1-alpha)*(1-p)*(1-h_e)*(1-uu-ud)^q) +
                            ud*alpha*(p/2 + (1-p)*h_p + (1-p)*(1-h_p)*(uu+ud)^q) +
                            du*(1-alpha)*(1 - (1-p)*(1-h_e)*(1-uu-ud)^q) - uu

                        F[2] = uu*alpha*p/2 +
                            ud*((1-alpha)*(1-p)*(h_e + (1-h_e)*(uu+ud)^q) +
                                alpha*(1-p/2) - alpha*(1-p)*(h_p + (1-h_p)*(uu+ud)^q)) +
                            dd*(1-alpha)*(1-p)*(h_e + (1-h_e)*(uu+ud)^q) - ud

                        F[3] = uu*(1-alpha)*(1-p)*(1-h_e)*(1-uu-ud)^q +
                            du*((1-alpha)*(1-p)*(1-h_e)*(1-uu-ud)^q +
                                alpha*(1 - p/2 - (1-p)*h_p - (1-p)*(1-h_p)*(1-uu-ud)^q)) +
                            dd*alpha*(p/2 + (1-p)*h_p) - du

                        F[4] = ud*(1-alpha)*(p + (1-p)*(1-h_e)*(1-(uu+ud)^q)) +
                            du*alpha*(p/2 + (1-p)*h_p + (1-p)*(1-h_p)*(1-uu-ud)^q) +
                            dd*(1 - alpha - (1-alpha)*(1-p)*(h_e + (1-h_e)*(uu+ud)^q) +
                            alpha*(1 - p/2 - (1-p)*h_p)) - dd
                    end

                    # Initial guess (choose reasonable starting point)
                    guess = [1., 0., 0., 0.]  #[0.5, 0.5, 0.5, 0.5]

                    sol = nlsolve(system!, guess; xtol=1e-15, ftol=1e-15)

                    if sol.f_converged
                        uu, ud, du, dd = sol.zero
                        cS = uu + ud
                        cSigma = uu + du
                        diss = ud + du

                        push!(cSUp, cS)
                        push!(cSigmaUp, cSigma)
                        push!(dissonance, diss)

                        # Hysteresis check (mimicking MATLAB length check)
                        # Here NLsolve returns only one root → need multi-start if you want multiple
                        # For now, keep single solution
                    end

                    # Write results immediately
                    open(export_file, "a") do io
                        @printf(io, "%.4f \t %.4f \t %.4f \t %.4f\n", p, cSUp[end], cSigmaUp[end], dissonance[end])
                    end
                end

                # Write hysteresis part if available
                for i in eachindex(cSup_h)
                    open(export_file, "a") do io
                        @printf(io, "%.4f \t %.4f \t %.4f \t %.4f\n", p_h[i], cSup_h[i], cSigmaUp_h[i], dissonance_h[i])
                    end
                end
            end
        end
    end
end