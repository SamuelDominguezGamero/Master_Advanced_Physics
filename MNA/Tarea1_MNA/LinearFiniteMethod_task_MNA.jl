# y'' = p(x)y' + q(x)y + r(x) -> ODE -> For this specific case: p(x) = 0; q(x) = -(1+x^2); r(x) = -1

# N --->  number of sub-divisions --> 2N + 1 POINTS ---> x_i = 0 +- i* l / N

# FINITE DIFFERENCE METHOD's EQUATION SYSTEM:
# (2+h**2 * q(x_i))* w_i - 1 * w_{i+1} = -h**2 * r(x_i) + (1)*w_{0}  -> for i = - (N-1) / N * l
# - w_{i-1} + (2+h^2 * q(x_i))* w_i - w_{i+1} = -h**2 * r(x_i) -> for a point within the inner grid (i = 2, 3, ..., N-1)
# - w_{i-1} + (2+h^2 * q(x_i))* w_i = -h**2 * r(x_i) + w_{N+1}  -> for i = - (N-1) / N * l

# Variable "i" should be renamed for simplicity in the code: i ---> [1, 2N + 1]

using LinearAlgebra

# y'' = p(x)y' + q(x)y + r(x)
p(x) = 0
q(x) = -(1.0 + x^2)
r(x) = -1.0

# Domain and Boundary Conditions
a, b = -1.0, 1.0
y_a, y_b = 0.0, 0.0

N = 100
M = 2 * N + 1

x = collect(range(a, b, length = M))
h = step(range(a, b, length = M))

x_int = x[2:end-1]    # Inner points of the Grid
N_int = length(x_int) # Number of inner points

dl = -1.0 .- (0.5 * h) .* p.(x_int[2:end])   # Subdiagonal
du = -1.0 .+ (0.5 * h) .* p.(x_int[1:end-1]) # Superdiagonal 
d  = 2.0 .+ (h^2) .* q.(x_int)               # Main diagonal

A = Tridiagonal(dl, d, du) # Tridiagonal Matrix

# Término general: -h^2 * r(x_i)
b_vec       = (-h^2) .* r.(x_int)
b_vec[1]   += (1.0 + 0.5 * h * p(x_int[1]))   * y_a
b_vec[end] += (1.0 - 0.5 * h * p(x_int[end])) * y_b

# Solving the System of Equations
w_int = A \ b_vec
y_sol = [y_a; w_int; y_b]


### SCRIPT DE PROTOTIPADO TERMINADO HASTA AQUÍ. Hemos resuelto la función del ejercicio.




















## Creación de una función modular de propósito general para Finite Difference Method:
using LinearAlgebra

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

    return x, y
end
