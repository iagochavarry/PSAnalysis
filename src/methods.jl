using Base: Float64
P_flux(B, V, θ, i, k) = V[i]*V[k]*B[i, k]*sin(θ[i]- θ[k])
Q_flux(B, V, θ, i, k) = -B[i,k]*V[i]^2 + V[i]*V[k]*B[i, k]*cos(θ[i]- θ[k])

P_eq(B, V, θ, i) = V[i]*(sum(V[k]*B[i, k]*sin(θ[i] - θ[k]) for k = 1:length(V)))
Q_eq(B, V, θ, i) = V[i]*(sum(-V[k]*B[i, k]*cos(θ[i] - θ[k]) for k = 1:length(V)))

P(B, V, θ) = [missing, P_eq(B, V, θ, 2), P_eq(B, V, θ, 3), P_eq(B, V, θ, 4), P_eq(B, V, θ, 5)]
Q(B, V, θ) = [missing, missing, Q_eq(B, V, θ, 3), Q_eq(B, V, θ, 4), Q_eq(B, V, θ, 5)]

J(B, V, θ) = [J_H(B, V, θ) J_N(B, V, θ); J_M(B, V, θ) J_L(B, V, θ)]

function J_H(B, V, θ)
    n_bus = length(B[:, 1])
    jac = Array{Any, 2}(undef, n_bus, n_bus)
    for i = 1:n_bus, k = 1:n_bus
        jac[i, k] = i==k ? J_H(B, V, θ, i) : J_H(B, V, θ, i, k)
    end
    return jac
end

J_H(B, V, θ, i, k) = - V[i]*V[k]*B[i, k]*cos(θ[i] - θ[k])
J_H(B, V, θ, i)    = V[i]*(sum(V[k]*B[i, k]*cos(θ[i] - θ[k]) for k in setdiff(collect(1:length(V)), [i])))

function J_N(B, V, θ)
    n_bus = length(B[:, 1])
    jac = Array{Any, 2}(undef, n_bus, n_bus)
    for i = 1:n_bus, k = 1:n_bus
        jac[i, k] = i==k ? J_N(B, V, θ, i) : J_N(B, V, θ, i, k)
    end
    return jac
end

J_N(B, V, θ, i, k) = V[i]*B[i, k]*sin(θ[i] - θ[k])
J_N(B, V, θ, i)    = sum(V[k]*B[i, k]*sin(θ[i] - θ[k]) for k in setdiff(collect(1:length(V)), [i]))

function J_M(B, V, θ)
    n_bus = length(B[:, 1])
    jac = Array{Any, 2}(undef, n_bus, n_bus)
    for i = 1:n_bus, k = 1:n_bus
        jac[i, k] = i==k ? J_M(B, V, θ, i) : J_M(B, V, θ, i, k)
    end
    return jac
end

J_M(B, V, θ, i, k) = -V[i]*V[k]*B[i, k]*sin(θ[i] - θ[k])
J_M(B, V, θ, i)    = V[i]*(sum(V[k]*B[i, k]*sin(θ[i] - θ[k]) for k in setdiff(collect(1:length(V)), [i])))

function J_L(B, V, θ)
    n_bus = length(B[:, 1])
    jac = Array{Any, 2}(undef, n_bus, n_bus)
    for i = 1:n_bus, k = 1:n_bus
        jac[i, k] = i==k ? J_L(B, V, θ, i) : J_L(B, V, θ, i, k)
    end
    return jac
end

J_L(B, V, θ, i, k) = - V[i]*B[i, k]*cos(θ[i] - θ[k])
J_L(B, V, θ, i)    = -2*V[i]*B[i,i] + sum(-V[k]*B[i, k]*cos(θ[i] - θ[k]) for k in setdiff(collect(1:length(V)), [i]))


function calc_jacobiano(B, V, θ)
    cols = [2,3,4,5,8,9,10]
    lines = [2,3,4,5,8,9,10]
    return J(B,V,θ)[lines, cols]
end

function has_converged(gap::Float64, tol::Float64)
    return gap > tol ? false : true 
end

function calc_P_flux(V::Vector, θ::Vector, B::Matrix)
    n_bus = size(B, 1)
    P = zeros(n_bus, n_bus)
    for i = 1:n_bus, k = 1:n_bus
        P[i, k] = i==k ? P_eq(B, V, θ, i) : P_flux(B, V, θ, i, k)
    end
    return P
end

function calc_Q_flux(V::Vector, θ::Vector, B::Matrix)
    n_bus = size(B, 1)
    Q = zeros(n_bus, n_bus)
    for i = 1:n_bus, k = 1:n_bus
        Q[i, k] = i==k ? Q_eq(B, V, θ, i) : Q_flux(B, V, θ, i, k)
    end
    return Q
end

function round_num(x::Number)
    return abs(x) <= 1e-7 ? 0.0 : round(x, digits = 7)
end

function round_num(x::Missing)
    return x
end

function newton_method(input::Input, 
                       θ_vars::Vector{Int}, V_vars::Vector{Int}, 
                       P_esp::Vector, Q_esp::Vector; 
                       iter_lim = 20, tol = 1e-5)
    YBUS = PSAnalysis.ybus(input)
    G = real.(YBUS)
    B = imag.(YBUS)
    ####
    n_bus = size(B, 1)
    var_cols =  vcat(θ_vars, n_bus .+ V_vars) # variables idxs in θ and V
    ####
    Δ  = Array{Any, 2}(missing, n_bus*2, iter_lim)
    gap = Vector{Float64}(undef, iter_lim)
    ####
    V = ones(5, iter_lim) # idx 1, 2 are fixed
    θ = zeros(5, iter_lim) # idx 1 is fixed
    iter = 1
    Δ[1:n_bus, iter] = P_esp .- P(B, V[:, iter], θ[:, iter])
    Δ[n_bus+1:end, iter] = Q_esp .- Q(B, V[:, iter], θ[:, iter])
    gap[iter] = maximum(skipmissing(abs.(Δ[:, iter])))
    converged = has_converged(gap[iter], tol)
    while (converged == false) && (iter+1 <= iter_lim)
        # Jacobiana
        jacob = calc_jacobiano(B, V[:, iter], θ[:, iter])
        # Calcular tamanho do passo
        δ = Float64.(jacob) \ Float64.(Δ[var_cols, iter])
        # Atualizar valores
        iter += 1
        θ[θ_vars, iter] = θ[θ_vars, iter-1] .+ δ[1:length(θ_vars)]
        V[V_vars, iter] = V[V_vars, iter-1] .+ δ[length(θ_vars)+ 1:end]
        Δ[1:n_bus, iter] = P_esp .- P(B, V[:, iter], θ[:, iter])
        Δ[n_bus+1:end, iter] = Q_esp .- Q(B, V[:, iter], θ[:, iter])
        gap[iter] = maximum(skipmissing(abs.(Δ[:, iter])))
        converged = has_converged(gap[iter], tol)
    end
    
    V_final = V[:, iter]
    θ_final = θ[:, iter]

    P_mat = calc_P_flux(V_final,θ_final, B)
    Q_mat = calc_Q_flux(V_final,θ_final, B)

    return  Results(
        round_num.(P_mat),
        round_num.(Q_mat),
        round_num.(sqrt.(P_mat.^2 .+ Q_mat.^2)),
        round_num.(V_final),
        round_num.(θ_final),
        AlgorithmResults(
            round_num.(V[:, 1:iter]),
            round_num.(θ[:, 1:iter]),
            round_num.(Δ[:, 1:iter]),
            round_num.(gap[1:iter])
        )
    )
end

function jacob_const_method(input::Input, 
                            θ_vars::Vector{Int}, V_vars::Vector{Int}, 
                            P_esp::Vector, Q_esp::Vector; 
                            iter_lim = 20, tol = 1e-5)
    YBUS = PSAnalysis.ybus(input)
    G = real.(YBUS)
    B = imag.(YBUS)
    ####
    n_bus = size(B, 1)
    var_cols =  vcat(θ_vars, n_bus .+ V_vars) # variables idxs in θ and V
    ####
    Δ  = Array{Any, 2}(undef, n_bus*2, iter_lim)
    gap = Vector{Float64}(undef, iter_lim)
    ####
    V = ones(5, iter_lim) # idx 1, 2 are fixed
    θ = zeros(5, iter_lim) # idx 1 is fixed
    iter = 1
    Δ[1:n_bus, iter] = P_esp .- P(B, V[:, iter], θ[:, iter])
    Δ[n_bus+1:end, iter] = Q_esp .- Q(B, V[:, iter], θ[:, iter])
    gap[iter] = maximum(skipmissing(abs.(Δ[:, iter])))
    converged = has_converged(gap[iter], tol)
    # Jacobiana
    jacob = calc_jacobiano(B, V[:, iter], θ[:, iter])
    while (converged == false) && (iter+1 <= iter_lim)
        # Calcular tamanho do passo
        δ = Float64.(jacob) \ Float64.(Δ[var_cols, iter])
        # Atualizar valores
        iter += 1
        θ[θ_vars, iter] = θ[θ_vars, iter-1] .+ δ[1:length(θ_vars)]
        V[V_vars, iter] = V[V_vars, iter-1] .+ δ[length(θ_vars)+ 1:end]
        Δ[1:n_bus, iter] = P_esp .- P(B, V[:, iter], θ[:, iter])
        Δ[n_bus+1:end, iter] = Q_esp .- Q(B, V[:, iter], θ[:, iter])
        gap[iter] = maximum(skipmissing(abs.(Δ[:, iter])))
        converged = has_converged(gap[iter], tol)
    end
    
    V_final = V[:, iter]
    θ_final = θ[:, iter]

    P_mat = calc_P_flux(V_final,θ_final, B)
    Q_mat = calc_Q_flux(V_final,θ_final, B)
    return  Results(
        round_num.(P_mat),
        round_num.(Q_mat),
        round_num.(sqrt.(P_mat.^2 .+ Q_mat.^2)),
        round_num.(V_final),
        round_num.(θ_final),
        AlgorithmResults(
            round_num.(V[:, 1:iter]),
            round_num.(θ[:, 1:iter]),
            round_num.(Δ[:, 1:iter]),
            round_num.(gap[1:iter])
        )
    )
end

function mix_method(input::Input, 
                            θ_vars::Vector{Int}, V_vars::Vector{Int}, 
                            P_esp::Vector, Q_esp::Vector; 
                            iter_lim = 20, tol = 1e-5, mix = 5e-1)
    YBUS = PSAnalysis.ybus(input)
    G = real.(YBUS)
    B = imag.(YBUS)
    ####
    n_bus = size(B, 1)
    var_cols =  vcat(θ_vars, n_bus .+ V_vars) # variables idxs in θ and V
    ####
    Δ  = Array{Any, 2}(undef, n_bus*2, iter_lim)
    gap = Vector{Float64}(undef, iter_lim)
    ####
    V = ones(5, iter_lim) # idx 1, 2 are fixed
    θ = zeros(5, iter_lim) # idx 1 is fixed
    iter = 1
    Δ[1:n_bus, iter] = P_esp .- P(B, V[:, iter], θ[:, iter])
    Δ[n_bus+1:end, iter] = Q_esp .- Q(B, V[:, iter], θ[:, iter])
    gap[iter] = maximum(skipmissing(abs.(Δ[:, iter])))
    converged = has_converged(gap[iter], tol)
    jacob = Float64[]
    # Jacobiana
    
    while (converged == false) && (iter+1 <= iter_lim)
        if (iter == 1) || (gap[iter] <= mix)
            @show iter, "entrei"
            jacob = calc_jacobiano(B, V[:, iter], θ[:, iter])
        end
        # Calcular tamanho do passo
        δ = Float64.(jacob) \ Float64.(Δ[var_cols, iter])
        # Atualizar valores
        iter += 1
        θ[θ_vars, iter] = θ[θ_vars, iter-1] .+ δ[1:length(θ_vars)]
        V[V_vars, iter] = V[V_vars, iter-1] .+ δ[length(θ_vars)+ 1:end]
        Δ[1:n_bus, iter] = P_esp .- P(B, V[:, iter], θ[:, iter])
        Δ[n_bus+1:end, iter] = Q_esp .- Q(B, V[:, iter], θ[:, iter])
        gap[iter] = maximum(skipmissing(abs.(Δ[:, iter])))
        converged = has_converged(gap[iter], tol)
    end
    
    V_final = V[:, iter]
    θ_final = θ[:, iter]

    P_mat = calc_P_flux(V_final,θ_final, B)
    Q_mat = calc_Q_flux(V_final,θ_final, B)
    return  Results(
        round_num.(P_mat),
        round_num.(Q_mat),
        round_num.(sqrt.(P_mat.^2 .+ Q_mat.^2)),
        round_num.(V_final),
        round_num.(θ_final),
        AlgorithmResults(
            round_num.(V[:, 1:iter]),
            round_num.(θ[:, 1:iter]),
            round_num.(Δ[:, 1:iter]),
            round_num.(gap[1:iter])
        )
    )
end

