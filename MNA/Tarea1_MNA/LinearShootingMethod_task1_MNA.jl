# ==============================================================================
# Sistema de EDOs: y'' = -(1+x^2)y - 1
# u = (y, v) donde v = y'
# Retorna (dy/dx, dv/dx)
# ==============================================================================
@inline function f(x::Float64, u::Tuple{Float64, Float64})
    y, v = u
    dydx = v # y' 
    dvdx = -(1.0 + x^2) * y - 1.0 # y''
    return (dydx, dvdx)
end

# ==============================================================================
# Rutina Runge-Kutta 4 (RK4) para un paso escalar/vectorial
# ==============================================================================
@inline function rk4_step(x::Float64, u::Tuple{Float64, Float64}, h::Float64)
    k1_y, k1_v = f(x, u)
    
    u_k2 = (u[1] + 0.5 * h * k1_y, u[2] + 0.5 * h * k1_v)
    k2_y, k2_v = f(x + 0.5 * h, u_k2)
    
    u_k3 = (u[1] + 0.5 * h * k2_y, u[2] + 0.5 * h * k2_v)
    k3_y, k3_v = f(x + 0.5 * h, u_k3)
    
    u_k4 = (u[1] + h * k3_y, u[2] + h * k3_v)
    k4_y, k4_v = f(x + h, u_k4)
    
    y_next = u[1] + (h / 6.0) * (k1_y + 2.0 * k2_y + 2.0 * k3_y + k4_y)
    v_next = u[2] + (h / 6.0) * (k1_v + 2.0 * k2_v + 2.0 * k3_v + k4_v)
    
    return (y_next, v_next)
end

# ==============================================================================
# Integrador global
# Discretización: x \in [-1, 1], con 2N subintervalos, 2N + 1 puntos
# ==============================================================================
function solve_rk4(N::Int, y0::Float64, v0::Float64)
    # Longitud total del dominio [-1, 1] es L = 2.0
    # Dividiendo en 2N intervalos simétricos respecto al origen:
    M = 2 * N
    x0 = -1.0
    xf = 1.0
    h = (xf - x0) / M

    # Vectores de almacenamiento
    x_vals = zeros(Float64, M + 1)
    y_vals = zeros(Float64, M + 1)
    v_vals = zeros(Float64, M + 1)

    # Condiciones iniciales
    x_vals[1] = x0
    y_vals[1] = y0
    v_vals[1] = v0

    u_current = (y0, v0)

    # Integración temporal explícita
    for i in 1:M
        x_current = x_vals[i]
        u_current = rk4_step(x_current, u_current, h)
        
        x_vals[i + 1] = x0 + i * h
        y_vals[i + 1] = u_current[1]
        v_vals[i + 1] = u_current[2]
    end

    return x_vals, y_vals, v_vals
end

# ==============================================================================
# Ejecución y comparación de ambos casos
# ==============================================================================

using Plots

# Resolver ambos casos
N = isempty(ARGS) ? 300 : parse(Int, ARGS[1])
x1, y1, _ = solve_rk4(N, 0.0, 0.0)
x2, y2, _ = solve_rk4(N, 0.0, 1.0)

y = y1 + y2*(-y1[end]/y2[end])

# Gráfica de la primera solución
plot(x1, y1, 
     label="Función auxiliar y_1", 
     lw=2, 
     color=:blue,
     xlabel="x", 
     ylabel="y(x)",
     legend=:topright,
     ls=:dash,
     grid=true)

# Superponer la segunda solución
plot!(x2, y2, 
      label="Función auxiliar y_2", 
      lw=2,
      color=:red, 
      ls=:dash)

plot!(x1, y, 
      label="Solución numérica al PVF", 
      lw=3, 
      color=:cyan)

# Guardar la figura en disco
savefig("soluciones_rk4_N=$N.pdf")
