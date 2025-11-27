using Symbolics
using NLsolve
using LinearAlgebra


function solvesystem(α_val, p_val, he_val, hp_val, q_val, initguess)
    function system!(F, x)
        uu, ud, du, dd = x
        
        F[1] = uu*(1 - α_val*p_val/2 - (1-α_val)*(1-p_val)*(1-he_val)*(1-uu-ud)^q_val) +
               ud*α_val*(p_val/2 + (1-p_val)*hp_val + (1-p_val)*(1-hp_val)*(uu+ud)^q_val) +
               du*(1-α_val)*(1 - (1-p_val)*(1-he_val)*(1-uu-ud)^q_val) - uu

        F[2] = uu*α_val*p_val/2 +
               ud*((1-α_val)*(1-p_val)*(he_val + (1-he_val)*(uu+ud)^q_val) + α_val*(1-p_val/2) - α_val*(1-p_val)*(hp_val + (1-hp_val)*(uu+ud)^q_val)) +
               dd*(1-α_val)*(1-p_val)*(he_val + (1-he_val)*(uu+ud)^q_val) - ud

        F[3] = uu*(1-α_val)*(1-p_val)*(1-he_val)*(1-uu-ud)^q_val +
               du*((1-α_val)*(1-p_val)*(1-he_val)*(1-uu-ud)^q_val + α_val*(1 - p_val/2 - (1-p_val)*hp_val - (1-p_val)*(1-hp_val)*(1-uu-ud)^q_val)) +
               dd*α_val*(p_val/2 + (1-p_val)*hp_val) - du

        F[4] = uu + ud + du + dd - 1
    end

    sol = nlsolve(system!, initguess)
    return sol.zero

end


function simplified_jacobian()
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


function checkstability(alpha_, p_, he_, hp_, q_, initialguess, J_simplified)
@variables uu ud du dd α p he hp q

    solution_point = solvesystem(alpha_, p_, he_, hp_, q_, initialguess)
    println("Fixed point solution: ", solution_point)

    # Substitute fixed point + parameters into symbolic Jacobian
    subs_dict = Dict(α => alpha_, p => p_, he => he_, hp => hp_, q => q_,
        uu => solution_point[1], ud => solution_point[2],
        du => solution_point[3], dd => solution_point[4])

    J_num = substitute(J_simplified, subs_dict) |> eval

    eigenvalues = eigvals(Matrix(J_num))
    println("Eigenvalues: ", eigenvalues)

    solution_point = solvesystem(alpha_, p_, he_, hp_, q_, initialguess)
    println("Fixed point solution: ", solution_point)

    # Calculate eigenvalues
    eigenvalues = eigvals(J_num)
    println("Eigenvalues: ", eigenvalues)

    # Stability analysis
    if all(abs.(eigenvalues) .< 1)
        println("Fixed point is STABLE")
        return true
    else
        println("Fixed point is UNSTABLE")
        return false
    end
    
end

jacobian_simplified = simplified_jacobian()

# checkstability(0.1, 0.244, 0.01, 0.0, 3, [0., 0., 0., 1.], jacobian_simplified)
# checkstability(0.1, 0.244, 0.01, 0.0, 3, [1, 0., 0., 0.], jacobian_simplified)
# # checkstability(0.1, 0.05, 0.0, 0.0, 3, [0.25, 0.25, 0.25, 0.25], jacobian_simplified)


# checkstability(0.1, 0.297, 0.01, 0.0, 3, [0., 0., 0., 1.], jacobian_simplified)
# checkstability(0.1, 0.297, 0.01, 0.0, 3, [1, 0., 0., 0.], jacobian_simplified)


checkstability(0.1, 0.2378, 0.01, 0.0, 3, [0., 0., 0., 1.], jacobian_simplified)
checkstability(0.1, 0.2378, 0.01, 0.0, 3, [1., 0., 0., 0.], jacobian_simplified)