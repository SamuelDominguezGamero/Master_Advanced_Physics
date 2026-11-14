using LinearAlgebra
using Plots
"""
    finite_difference_bvp(p, q, r, a, b, y_a, y_b, N)

Resuelve el PVF lineal 1D de segundo orden:
    y'' = p(x)y' + q(x)y + r(x),   x ∈ [a, b]
    y(a) = y_a,  y(b) = y_b

Parámetros:
- `p, q, r`: Funciones univariables f(x).
- `a, b`: Límites del dominio espacial.
- `y_a, y_b`: Condiciones de contorno tipo Dirichlet.
- `N`: Número de subdivisiones (genera 2N + 1 puntos totales en la malla).

Devuelve:
- `x`: Vector con la malla completa de puntos (longitud 2N + 1).
- `y`: Vector con la solución numérica evaluada en `x`.
"""
function finite_difference_bvp(p::F1, q::F2, r::F3, 
                               a::Real, b::Real, 
                               y_a::Real, y_b::Real, 
                               N::Integer) where {F1, F2, F3}
    M = 2 * N + 1
    r_grid = range(a, b, length = M)
    h = step(r_grid)
    x = collect(r_grid)

    # Nodos interiores
    x_int = x[2:end-1]
    N_int = length(x_int)

    # Construcción de las diagonales
    dl = -1.0 .- (0.5 * h) .* p.(x_int[2:end])
    du = -1.0 .+ (0.5 * h) .* p.(x_int[1:end-1])
    d  =  2.0 .+ (h^2)      .* q.(x_int)

    A = Tridiagonal(dl, d, du)

    # Vector de términos independientes
    rhs = (-h^2) .* r.(x_int)
    rhs[1]   += (1.0 + 0.5 * h * p(x_int[1]))   * y_a
    rhs[end] += (1.0 - 0.5 * h * p(x_int[end])) * y_b

    # Resolución O(N) y ensamblado del perfil completo
    w_int = A \ rhs
    y = [y_a; w_int; y_b]

    return x, y, h
end


a = -1.0
b =  1.0
y_a = 0.0
y_b = 0.0
N = 100

p(x) = 0.0
q(x) = -(1.0 + x^2)
r(x) = -1.0


#  Solución de referencia (malla muy fina)
N_ref = 5000
x_ref, y_ref, _ = finite_difference_bvp(p, q, r, a, b, y_a, y_b, N_ref)

# Barrido sobre diferentes resoluciones
N_vals = [5, 10, 20, 50, 200]
h_vals = Float64[]
err_L2 = Float64[]
err_Linf = Float64[]

# Inicialización de la figura 1
plt_profiles = plot(
    title = "Comparativa visual de la solución numérica y(x)",
    xlabel = "x",
    ylabel = "y(x)",
    legend = :topright,
    grid = true
)

for N in N_vals
    x_n, y_n, h_n = finite_difference_bvp(p, q, r, a, b, y_a, y_b, N)
    
    # Añadir traza al gráfico de perfiles
    plot!(plt_profiles,
    	  x_n, y_n, 
    	  label = "N = $N (h = $(round(h_n, digits=4)))", 
    	  lw = 1.5)
    	   
    # Comparar con la solución de referencia en los puntos comunes
    # Se muestrea la malla de referencia en los índices exactos correspondientes
    stride = div(N_ref, N)
    y_ref_sampled = y_ref[1:stride:end]
    
    push!(h_vals, h_n) # guarda h_n en h_vals
    push!(err_Linf, norm(y_n .- y_ref_sampled, Inf)) # De todos los errores sampleados, guarda el más grande en valor absoluto 
end

# Traza de la solución de referencia
plot!(plt_profiles, x_ref, y_ref, label = "Referencia (N = $N_ref)", color = :black, linestyle = :dash, lw = 2)

savefig(plt_profiles, "perfiles_convergencia.pdf")


# -------------------------------------------------------------
# Figura 2: Curva de orden de convergencia (log-log)
# -------------------------------------------------------------
plt_error = plot(
    h_vals, err_Linf,
    xscale = :log10,
    yscale = :log10,
    marker = :circle,
    lw = 2,
    label = "Error ||y_N - y_ref||_inf",
    xlabel = "Paso de malla h",
    ylabel = "Error relativo a referencia",
    title = "Orden de convergencia (Diferencias Centradas)",
    legend = :bottomright,
    grid = true
)

# Línea de referencia teórica O(h^2)
c_fit = err_Linf[end] / (h_vals[end]^2)
plot!(plt_error, h_vals, c_fit .* (h_vals .^ 2), 
      label = "Pendiente teórica O(h²)", 
      linestyle = :dot, 
      color = :black, 
      lw = 1.5)

savefig(plt_error, "orden_convergencia_loglog.pdf")
