using Pkg
Pkg.add("Parameters")
Pkg.add("BenchmarkTools")
using BenchmarkTools
include("src/PSAnalysis.jl")

# Barras
buses = [1, 2, 3, 4, 5]
# Linhas
lines = [1, 2, 3, 3, 5]
# De qual barra sai cada linha
from = [1,3,4,3,4]
# Para qual barra vai cada linha
to = [3,4,2,5,5]
# Resistencia da linha
r_line = [0.0, 0.0, 0.0, 0.0, 0.0]
# Reatância da linha
x_line = [0.02, 0.04, 0.02, 0.04, 0.04]
# Capactor Shunt em cada Barra
capacitor_bank = Complex{Float64}[0,0,0,0, 4.4im]

# Definindo Input do sistema
input = PSAnalysis.Input(
    buses, lines, from, to, 
    r_line, x_line, capacitor_bank = capacitor_bank,
    type_line = "impedance"
)

PSAnalysis.ybus(input)

# Indíce das variáveis de ângulo
θ_vars = [2,3,4,5]
# Indíce das variáveis de tensão
V_vars = [3,4,5]
# Potência Ativa Especificada
P_esp = [missing, 5.5,0,0,-8]
# Potência Reativa Especificada
Q_esp = [missing, missing, 0,0,-3]

@time results_nr = PSAnalysis.newton_method(input, θ_vars, V_vars, P_esp, Q_esp, tol = 1e-5);

@time results_jc = PSAnalysis.jacob_const_method(input, θ_vars, V_vars, P_esp, Q_esp, tol = 1e-5)

@time results_mx = PSAnalysis.mix_method(input, θ_vars, V_vars, P_esp, Q_esp, tol = 1e-5, mix = 0.6);


#### Realizando Benchmark ####
@benchmark bench_nr = PSAnalysis.newton_method(input, θ_vars, V_vars, P_esp, Q_esp, tol = 1e-5)

@benchmark bench_jc = PSAnalysis.jacob_const_method(input, θ_vars, V_vars, P_esp, Q_esp, tol = 1e-5)

@benchmark bench_mx = PSAnalysis.mix_method(input, θ_vars, V_vars, P_esp, Q_esp, tol = 1e-5, mix = 0.6)


# Resultado final do método de newton
results_nr.P
results_nr.Q
results_nr.S
results_nr.V
results_nr.θ
## Valores durante as iterações do método de newton
results_nr.algorithm_results.V
results_nr.algorithm_results.θ
results_nr.algorithm_results.Δ
results_nr.algorithm_results.gap


# Resultado final do método do jacobiano constante
results_jc.P
results_jc.Q
results_jc.S
results_jc.V
results_jc.θ
## Valores durante as iterações do jacobiano constante
results_jc.algorithm_results.V
results_jc.algorithm_results.θ
results_jc.algorithm_results.Δ
results_jc.algorithm_results.gap


############## Plot Iterações NEWTON ##############
p0 = plot(size = (1200, 700), title = "Newton-Rapson Iterações", axis = false, grid = false, border = :none)
p1 = plot(size = (1200, 700), title = "Tensão [p.u.]")
plot!(results_nr.algorithm_results.V[3,:], label = "V3", color=:gray70)
plot!(results_nr.algorithm_results.V[4,:], label = "V4", color=:gray50)
plot!(results_nr.algorithm_results.V[5,:], label = "V5", color=:gray30)

p2 = plot(size = (1200, 700), title = "Ângulo [graus]", xaxis = "Iterações")
plot!(results_nr.algorithm_results.θ[2,:], label = "θ2", color =:gray70)
plot!(results_nr.algorithm_results.θ[3,:], label = "θ3", color =:gray60)
plot!(results_nr.algorithm_results.θ[4,:], label = "θ4", color =:gray50)
plot!(results_nr.algorithm_results.θ[5,:], label = "θ5", color =:gray40)

p3 = plot(size = (1200, 700), title = "Gap", xaxis = "Iterações")
plot!(results_nr.algorithm_results.gap, label = "", color =:gray40)

pl = plot(p0, p1, p2, p3, layout = @layout([a{0.01h};[b c]; d]))



############## Plot Iterações JACOBIANO CONSTANTE ##############

results_jc.algorithm_results.V

p0 = plot(size = (1200, 700), title = "Jacobiano Constante Iterações", axis = false, grid = false, border = :none)
p1 = plot(size = (1200, 700), title = "Tensão [p.u.]")
plot!(results_jc.algorithm_results.V[3,:], label = "V3", color=:gray70)
plot!(results_jc.algorithm_results.V[4,:], label = "V4", color=:gray50)
plot!(results_jc.algorithm_results.V[5,:], label = "V5", color=:gray30)

p2 = plot(size = (1200, 700), title = "Ângulo [graus]", xaxis = "Iterações")
plot!(results_jc.algorithm_results.θ[2,:], label = "θ2", color =:gray70)
plot!(results_jc.algorithm_results.θ[3,:], label = "θ3", color =:gray60)
plot!(results_jc.algorithm_results.θ[4,:], label = "θ4", color =:gray50)
plot!(results_jc.algorithm_results.θ[5,:], label = "θ5", color =:gray40)

p3 = plot(size = (1200, 700), title = "Gap", xaxis = "Iterações")
plot!(results_jc.algorithm_results.gap, label = "", color =:gray40)

pl = plot(p0, p1, p2, p3, layout = @layout([a{0.01h};[b c]; d]))

