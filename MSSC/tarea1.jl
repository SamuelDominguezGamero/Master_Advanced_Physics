using Random
using SpecialFunctions
using Statistics

"""
    estimar_pi_mc(N::Int, L::Float64, reps::Int; verbose::Bool=false)

Calcula estimaciones de Pi en 5D mediante Monte Carlo multihilo.
Devuelve: `(buffer_pi, media_pi, std_pi)`
"""
function estimar_pi_mc(N::Int, L::Float64, reps::Int; verbose::Bool=false)
    factorial_2_5 = gamma(3.5)
    r2_max = L * L
    dos_L = 2.0 * L
    
    # Buffer prealocado thread-safe
    buffer_pi = zeros(Float64, reps)
    
    # Distribución paralela sobre las realizaciones (semillas 1:reps)
    Threads.@threads for i in 1:reps
        rng = Xoshiro(i)
        puntos_dentro = 0
        
        for _ in 1:N
            r2 = 0.0
            @inline for _ in 1:5
                val = (rand(rng) * dos_L) - L
                r2 += val * val
            end
            if r2 < r2_max
                puntos_dentro += 1
            end
        end
        
        ratio = puntos_dentro / N
        buffer_pi[i] = ((32.0 * ratio) * factorial_2_5)^(0.4)
    end
    
    media_pi = mean(buffer_pi)
    std_pi   = std(buffer_pi)
    
    if verbose
        println("Buffer de Pi: ", buffer_pi)
        println("Media: ", media_pi, " ± ", std_pi)
    end
    
    return buffer_pi, media_pi, std_pi
end



# Primera Ejecución, para el ejercicio 1a
println("Primera iteración del experimento") 
buffer_pi_1, media_pi_1, std_pi_1 = estimar_pi_mc(1_000, 1.0, 10; verbose=true)



println("------------------------------------")


# Barrido logarítmico en N
println("Repetición del experimento para varios valores de N")
valores_N = [20^3, 10^4, 20^4, 10^5, 10^6, 10^7, 10^8]
L = 1.0
reps = 50

desviaciones = Float64[]

for N in valores_N
    _, _, sigma = estimar_pi_mc(N, L, reps; verbose=false)
    push!(desviaciones, sigma)
    println("N = $(rpad(N, 8)) -> std = $sigma")
end

# Comprobación de ley de escala (Monte Carlo estándar: sigma ~ 1 / sqrt(N))
# En escala log-log, la pendiente debe aproximarse a -0.5:
pendientes = diff(log.(desviaciones)) ./ diff(log.(valores_N))
println("Pendiente log-log estimada: ", mean(pendientes))


# -----------------------------------------------------------------

using Plots

# 1. Ajuste lineal en escala logarítmica: log(sigma) = a + b * log(N)
log_N     = log.(valores_N)
log_sigma = log.(desviaciones)

# Matriz de diseño para mínimos cuadrados [1 log(N)]
X = hcat(ones(length(log_N)), log_N)
beta = X \ log_sigma
ordenada = beta[1]
pendiente = beta[2]

# Recta de ajuste evaluada en los puntos
sigma_ajuste = exp.(ordenada .+ pendiente .* log_N)

# 2. Generación de la gráfica log-log
plot(
    valores_N, 
    desviaciones, 
    seriestype = :scatter, 
    xscale = :log10, 
    yscale = :log10,
    label = "Datos MC (reps = $reps)", 
    marker = (:circle, 6, :blue),
    xlabel = "Número de muestras (N)", 
    ylabel = "Desviación estándar σ(N)",
    title = "Convergencia Monte Carlo: σ(N) ∝ N^($round(pendiente, digits=3))",
    legend = :topright,
    grid = true,
    minorgrid = true
)

# Línea del ajuste empírico
plot!(
    valores_N, 
    sigma_ajuste, 
    line = (:dash, 2, :red), 
    label = "Ajuste lineal (pendiente = $(round(pendiente, digits=3)))"
)

# Curva de referencia teórica N^(-1/2) anclada al primer punto
sigma_teorica = desviaciones[1] .* sqrt.(valores_N[1] ./ valores_N)
plot!(
    valores_N, 
    sigma_teorica, 
    line = (:dot, 1.5, :gray), 
    label = "Referencia teórica N^(-1/2)"
)

# Guardar figura
savefig("convergencia_mc_pi_log.png")

